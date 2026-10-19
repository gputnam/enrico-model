#include "orchestration/Generator.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <map>
#include <mutex>
#include <sstream>
#include <thread>

#include "core/Config.h"
#include "core/PDG.h"
#include "core/Parallel.h"
#include "orchestration/ModelFactory.h"
#include "orchestration/XSecTableBuilder.h"

namespace nuflow {

Generator::Generator(const YAML::Node& config) : config_(YAML::Clone(config)) {
  const YAML::Node run = config_["run"];
  seed_ = Get<uint64_t>(run, "seed", 1);
  threads_ = Get<int>(run, "threads", static_cast<int>(std::max(1u, std::thread::hardware_concurrency())));
  batchSize_ = Get<int>(run, "batch_size", 5000);
  safetyFactor_ = Get<double>(run, "safety_factor", 1.2);

  nuPdg_ = config_["beam"]["pdg"].as<int>();
  if (!pdg::IsNeutrino(nuPdg_)) throw std::runtime_error("beam.pdg must be a neutrino PDG code");
  flux_ = MakeFlux(config_["beam"]["flux"]);
  target_ = MakeTarget(config_["target"]);
  channels_ = MakeChannels(config_["channels"]);
  fsi_ = MakeFSI(config_["fsi"]);

  // Stage-1 processes: every (channel, nucleus) pair that can occur.
  std::vector<Process> processes;
  for (size_t ic = 0; ic < channels_.size(); ++ic) {
    const YAML::Node chNode = config_["channels"][ic];
    for (const auto& comp : target_.components) {
      if (!channels_[ic]->Applies(nuPdg_, comp.nucleus)) continue;
      Process p;
      p.channel = channels_[ic].get();
      p.nucleus = comp.nucleus;
      p.count = comp.count;
      p.multiplicity = channels_[ic]->Multiplicity(nuPdg_, comp.nucleus);
      YAML::Emitter em;
      em << YAML::Flow << chNode;
      p.physicsKey = em.c_str();
      processes.push_back(p);
    }
  }
  if (processes.empty()) throw std::runtime_error("no channel applies to this neutrino and target");

  const YAML::Node xs = config_["xsec_table"];
  XSecTableOptions opt;
  const double emin = std::max(flux_->EMin(), Get<double>(xs, "emin", 0.05));
  const double emax = flux_->EMax();
  opt.energies = flux_->IsMonoEnergetic() ? std::vector<double>{flux_->EMin()}
                                          : MakeEnergyGrid(emin, emax, Get<int>(xs, "nodes", 150));
  opt.samplesPerNode = Get<long>(xs, "samples", 200000);
  opt.seed = Get<uint64_t>(xs, "seed", 12345);
  opt.threads = threads_;
  const std::string cache = Get<std::string>(xs, "cache", "data/cache");
  if (cache != "none") {
    // Relative cache directories live under $NUFLOW_DIR (or the source tree), not the cwd.
    const char* dir = std::getenv("NUFLOW_DIR");
    opt.cacheDir = (cache[0] == '/') ? cache : std::string(dir ? dir : NUFLOW_SOURCE_DIR) + "/" + cache;
  }
  BuildXSecTables(nuPdg_, processes, opt);

  selector_ = std::make_unique<ChannelSelector>(std::move(processes));
  sampler_ = std::make_unique<EnergySampler>(*flux_, *selector_);
}

Generator::~Generator() = default;

RunInfo Generator::Info() const {
  RunInfo info;
  YAML::Emitter em;
  em << config_;
  info.configYaml = em.c_str();
  info.nuPdg = nuPdg_;
  info.target = target_;
  info.flux = flux_.get();
  info.selector = selector_.get();
  info.fsiModel = fsi_->Name();
  info.seed = seed_;
  info.fluxAveragedSigma = sampler_->FluxAveragedSigma();
  if (const auto* h = dynamic_cast<const HistogramFlux*>(flux_.get())) {
    info.fluxIntegral = h->Integral();
    info.fluxUnits = h->Units();
  }
  return info;
}

long Generator::GenerateOne(Random& rng, Event& ev, long long* nOverweight, double* maxRatio) const {
  const double enu = sampler_->Sample(rng);
  const size_t ip = selector_->Select(enu, rng);
  const Process& proc = selector_->Get(ip);
  const double envelope = safetyFactor_ * proc.table.WMax(enu);
  long trials = 0;
  for (;;) {
    ++trials;
    const double w = proc.channel->Trial(nuPdg_, enu, proc.nucleus, ev, rng);
    if (w <= 0) continue;
    if (maxRatio) *maxRatio = std::max(*maxRatio, w / envelope);
    if (w > envelope) {  // envelope violated: keep the event (least-biased choice) but count it
      if (nOverweight) ++*nOverweight;
      break;
    }
    if (rng.Uniform() * envelope < w) break;
  }
  ev.channel = proc.channel->ProcessId();
  ev.weight = 1.0;
  ev.totalXSec = selector_->TotalSigma(enu);
  ev.processXSec = proc.table.Sigma(enu);
  fsi_->Propagate(proc.nucleus, ev, rng);
  return trials;
}

RunSummary Generator::Run(long long nEvents, const std::vector<OutputWriter*>& writers) {
  const auto t0 = std::chrono::steady_clock::now();
  const RunInfo info = Info();
  for (auto* w : writers) w->Begin(info);

  const long long nBatches = (nEvents + batchSize_ - 1) / batchSize_;
  const long long maxInFlight = 4LL * threads_ + 4;
  std::mutex mtx;
  std::condition_variable cv;
  std::map<long long, std::vector<Event>> ready;
  long long nextToWrite = 0;
  bool failed = false;

  std::atomic<long long> totalTrials{0}, totalOver{0};
  double maxRatio = 0;

  // Single writer thread: consumes batches strictly in order.
  std::thread writerThread([&] {
    long long written = 0;
    while (nextToWrite < nBatches) {
      std::vector<Event> batch;
      {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&] { return failed || ready.count(nextToWrite); });
        if (failed) return;
        batch = std::move(ready[nextToWrite]);
        ready.erase(nextToWrite);
      }
      for (const auto& ev : batch)
        for (auto* w : writers) w->Write(ev);
      written += static_cast<long long>(batch.size());
      {
        std::lock_guard<std::mutex> lock(mtx);
        ++nextToWrite;
      }
      cv.notify_all();
      if (nBatches >= 20 && nextToWrite % (nBatches / 10) == 0)
        std::printf("[nuflow] %lld / %lld events\n", written, nEvents), std::fflush(stdout);
    }
  });

  try {
    ParallelFor(static_cast<size_t>(nBatches), threads_, [&](size_t b, int) {
      {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&] { return failed || static_cast<long long>(b) < nextToWrite + maxInFlight; });
        if (failed) return;
      }
      const long long first = static_cast<long long>(b) * batchSize_;
      const long long n = std::min<long long>(batchSize_, nEvents - first);
      std::vector<Event> batch(n);
      Random rng(seed_, b);
      long long trials = 0, over = 0;
      double ratio = 0;
      for (long long i = 0; i < n; ++i) {
        batch[i].particles.reserve(12);
        trials += GenerateOne(rng, batch[i], &over, &ratio);
        batch[i].number = first + i;
      }
      totalTrials += trials;
      totalOver += over;
      {
        std::lock_guard<std::mutex> lock(mtx);
        maxRatio = std::max(maxRatio, ratio);
        ready.emplace(static_cast<long long>(b), std::move(batch));
      }
      cv.notify_all();
    });
  } catch (...) {
    {
      std::lock_guard<std::mutex> lock(mtx);
      failed = true;
    }
    cv.notify_all();
    writerThread.join();
    throw;
  }
  writerThread.join();

  RunSummary s;
  s.nEvents = nEvents;
  s.nTrials = totalTrials;
  s.nOverweight = totalOver;
  s.maxOverweightRatio = maxRatio;
  s.wallSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  for (auto* w : writers) w->End(s);
  return s;
}

}  // namespace nuflow
