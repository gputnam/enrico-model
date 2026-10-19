#include "orchestration/XSecTableBuilder.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "core/Config.h"
#include "core/Parallel.h"
#include "core/Random.h"

namespace nuflow {

namespace {

constexpr uint64_t kCacheMagic = 0x6e75666c6f777831ULL;  // "nuflowx1"
// Bump whenever physics code changes the trial weights for an unchanged configuration, so that
// stale cached tables are never reused.
constexpr int kCacheVersion = 2;

std::string CacheKey(int nuPdg, const Process& p, const XSecTableOptions& opt) {
  std::ostringstream s;
  s.precision(17);
  s << "v" << kCacheVersion << "|nu" << nuPdg << "|" << p.channel->Describe() << "|" << p.physicsKey << "|A" << p.nucleus.pdg << "|m"
    << p.nucleus.mass << "|kF" << p.nucleus.kFermiProton << "," << p.nucleus.kFermiNeutron << "|Eb"
    << p.nucleus.removalEnergy << "|N" << opt.samplesPerNode << "|seed" << opt.seed << "|E";
  for (double e : opt.energies) s << e << ",";
  return s.str();
}

bool ReadCache(const std::string& path, uint64_t hash, XSecTable& t) {
  std::ifstream in(path, std::ios::binary);
  if (!in) return false;
  uint64_t magic = 0, h = 0, n = 0;
  in.read(reinterpret_cast<char*>(&magic), 8);
  in.read(reinterpret_cast<char*>(&h), 8);
  in.read(reinterpret_cast<char*>(&n), 8);
  if (!in || magic != kCacheMagic || h != hash || n > 100000) return false;
  for (auto* v : {&t.energy, &t.sigma, &t.sigmaErr, &t.wMax}) {
    v->resize(n);
    in.read(reinterpret_cast<char*>(v->data()), static_cast<std::streamsize>(n * sizeof(double)));
  }
  return static_cast<bool>(in);
}

void WriteCache(const std::string& path, uint64_t hash, const XSecTable& t) {
  const std::string tmp = path + ".tmp";
  {
    std::ofstream out(tmp, std::ios::binary);
    const uint64_t n = t.energy.size();
    out.write(reinterpret_cast<const char*>(&kCacheMagic), 8);
    out.write(reinterpret_cast<const char*>(&hash), 8);
    out.write(reinterpret_cast<const char*>(&n), 8);
    for (const auto* v : {&t.energy, &t.sigma, &t.sigmaErr, &t.wMax})
      out.write(reinterpret_cast<const char*>(v->data()), static_cast<std::streamsize>(n * sizeof(double)));
  }
  std::filesystem::rename(tmp, path);
}

}  // namespace

std::vector<double> MakeEnergyGrid(double emin, double emax, int n) {
  if (emax <= emin || n < 2) return {emax};
  std::vector<double> g(n);
  const double r = std::log(emax / emin);
  for (int i = 0; i < n; ++i) g[i] = emin * std::exp(r * i / (n - 1));
  g.back() = emax;
  return g;
}

void BuildXSecTables(int nuPdg, std::vector<Process>& processes, const XSecTableOptions& opt) {
  const size_t nE = opt.energies.size();
  for (size_t ip = 0; ip < processes.size(); ++ip) {
    Process& proc = processes[ip];
    const std::string key = CacheKey(nuPdg, proc, opt);
    const uint64_t hash = HashString(key);
    std::string cachePath;
    if (!opt.cacheDir.empty()) {
      std::filesystem::create_directories(opt.cacheDir);
      char name[64];
      std::snprintf(name, sizeof(name), "xsec_%016llx.bin", static_cast<unsigned long long>(hash));
      cachePath = (std::filesystem::path(opt.cacheDir) / name).string();
      if (ReadCache(cachePath, hash, proc.table)) {
        if (opt.verbose)
          std::printf("[nuflow] %s on %s: cross-section table loaded from cache %s\n", proc.channel->Name().c_str(),
                      proc.nucleus.name.c_str(), cachePath.c_str());
        continue;
      }
    }

    const auto t0 = std::chrono::steady_clock::now();
    XSecTable& t = proc.table;
    t.energy = opt.energies;
    t.sigma.assign(nE, 0);
    t.sigmaErr.assign(nE, 0);
    t.wMax.assign(nE, 0);
    ParallelFor(nE, opt.threads, [&](size_t ie, int) {
      Random rng(opt.seed ^ hash, 0x5EED0000ULL + ie);
      Event scratch;
      scratch.particles.reserve(16);
      double sum = 0, sum2 = 0, wmax = 0;
      for (long i = 0; i < opt.samplesPerNode; ++i) {
        const double w = proc.channel->Trial(nuPdg, t.energy[ie], proc.nucleus, scratch, rng);
        sum += w;
        sum2 += w * w;
        wmax = std::max(wmax, w);
      }
      const double n = static_cast<double>(opt.samplesPerNode);
      const double mean = sum / n;
      const double var = std::max(0.0, sum2 / n - mean * mean);
      t.sigma[ie] = proc.multiplicity * mean;
      t.sigmaErr[ie] = proc.multiplicity * std::sqrt(var / n);
      t.wMax[ie] = wmax;
    });
    if (opt.verbose) {
      const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      std::printf("[nuflow] %s on %s: built %zu-node cross-section table (%ld samples/node) in %.1f s\n",
                  proc.channel->Name().c_str(), proc.nucleus.name.c_str(), nE, opt.samplesPerNode, dt);
    }
    if (!cachePath.empty()) WriteCache(cachePath, hash, t);
  }
}

}  // namespace nuflow
