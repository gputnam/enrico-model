// nuflow-gen: generate neutrino interaction events.
//
//   nuflow-gen --config config/minerva_le_fhc_numu_CH.yaml [-n N] [-o out.root] [--seed S] [--threads T]
//              [--mono E] [--nu PDG] [--target "1000060120:1,1000010010:1"] [--hepmc out.hepmc3]

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>

#include "core/Config.h"
#include "orchestration/Generator.h"
#include "orchestration/RootTreeWriter.h"
#ifdef NUFLOW_HAS_HEPMC3
#include "orchestration/NuHepMCWriter.h"
#endif

using namespace nuflow;

static void Usage() {
  std::printf(
      "usage: nuflow-gen --config FILE [options]\n"
      "  -n N              number of events (default: run.events)\n"
      "  -o FILE           ROOT output file (default: run.output)\n"
      "  --hepmc FILE      also write NuHepMC (HepMC3 ASCII) output\n"
      "  --seed S          random seed\n"
      "  --threads T       worker threads\n"
      "  --compression C   ROOT compression setting (default 505 = ZSTD level 5)\n"
      "  --mono E          mono-energetic beam of energy E [GeV] (overrides the flux)\n"
      "  --nu PDG          neutrino species\n"
      "  --target SPEC     target mixture, e.g. 1000060120:1,1000010010:1\n");
}

int main(int argc, char** argv) {
  std::string configPath, hepmcPath;
  YAML::Node overrides(YAML::NodeType::Map);
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) {
        std::fprintf(stderr, "missing value for %s\n", a.c_str());
        std::exit(2);
      }
      return argv[++i];
    };
    if (a == "--config" || a == "-c") configPath = next();
    else if (a == "-n") overrides["run"]["events"] = std::stoll(next());
    else if (a == "-o") overrides["run"]["output"] = next();
    else if (a == "--hepmc") hepmcPath = next();
    else if (a == "--seed") overrides["run"]["seed"] = std::stoull(next());
    else if (a == "--threads") overrides["run"]["threads"] = std::stoi(next());
    else if (a == "--compression") overrides["run"]["compression"] = std::stoi(next());
    else if (a == "--mono") overrides["mono"] = std::stod(next());
    else if (a == "--nu") overrides["nu"] = std::stoi(next());
    else if (a == "--target") overrides["target"] = next();
    else if (a == "-h" || a == "--help") return Usage(), 0;
    else {
      std::fprintf(stderr, "unknown argument %s\n", a.c_str());
      return Usage(), 2;
    }
  }
  if (configPath.empty()) return Usage(), 2;

  try {
    YAML::Node config = LoadConfig(configPath);
    for (const auto& kv : overrides["run"]) config["run"][kv.first.as<std::string>()] = kv.second;
    if (overrides["mono"]) {
      YAML::Node f;
      f["type"] = "mono";
      f["energy"] = overrides["mono"].as<double>();
      config["beam"]["flux"] = f;
    }
    if (overrides["nu"]) config["beam"]["pdg"] = overrides["nu"].as<int>();
    if (overrides["target"]) config["target"] = overrides["target"].as<std::string>();

    const long long nEvents = Get<long long>(config["run"], "events", 10000);
    const std::string output = Get<std::string>(config["run"], "output", "nuflow.root");

    Generator gen(config);
    std::printf("[nuflow] flux-averaged total cross section: %.5g x 1e-38 cm^2 per target unit\n",
                gen.Sampler().FluxAveragedSigma());

    std::vector<std::unique_ptr<OutputWriter>> writers;
    writers.push_back(std::make_unique<RootTreeWriter>(output, Get<int>(config["run"], "compression", 505)));
    if (!hepmcPath.empty()) {
#ifdef NUFLOW_HAS_HEPMC3
      writers.push_back(std::make_unique<NuHepMCWriter>(hepmcPath));
#else
      throw std::runtime_error("built without HepMC3 support");
#endif
    }
    std::vector<OutputWriter*> raw;
    for (auto& w : writers) raw.push_back(w.get());

    const RunSummary s = gen.Run(nEvents, raw);
    std::printf("[nuflow] generated %lld events in %.2f s (%.0f events/s), %lld trials (efficiency %.3f)\n",
                s.nEvents, s.wallSeconds, s.nEvents / s.wallSeconds, s.nTrials,
                double(s.nEvents) / double(s.nTrials));
    std::printf("[nuflow] envelope violations: %lld (%.2e of events); max w/w_envelope over all trials = %.3f\n", s.nOverweight,
                double(s.nOverweight) / double(s.nEvents), s.maxOverweightRatio);
    std::printf("[nuflow] output written to %s\n", output.c_str());
  } catch (const std::exception& e) {
    std::fprintf(stderr, "nuflow-gen: error: %s\n", e.what());
    return 1;
  }
  return 0;
}
