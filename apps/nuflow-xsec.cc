// nuflow-xsec: print the tabulated total cross sections for a configuration.
//
//   nuflow-xsec --config FILE [--nu PDG] [--target SPEC] [--emin E --emax E --nodes N]

#include <cstdio>
#include <string>

#include "core/Config.h"
#include "orchestration/Generator.h"

using namespace nuflow;

int main(int argc, char** argv) {
  std::string configPath;
  YAML::Node over(YAML::NodeType::Map);
  for (int i = 1; i + 1 < argc; i += 2) {
    const std::string a = argv[i], v = argv[i + 1];
    if (a == "--config") configPath = v;
    else if (a == "--nu") over["nu"] = std::stoi(v);
    else if (a == "--target") over["target"] = v;
    else if (a == "--emin") over["emin"] = std::stod(v);
    else if (a == "--emax") over["emax"] = std::stod(v);
    else if (a == "--nodes") over["nodes"] = std::stoi(v);
    else if (a == "--samples") over["samples"] = std::stol(v);
  }
  if (configPath.empty()) {
    std::fprintf(stderr, "usage: nuflow-xsec --config FILE [--nu PDG] [--target SPEC] [--emin E] [--emax E] "
                         "[--nodes N] [--samples N]\n");
    return 2;
  }
  try {
    YAML::Node config = LoadConfig(configPath);
    if (over["nu"]) config["beam"]["pdg"] = over["nu"];
    if (over["target"]) config["target"] = over["target"];
    if (over["emin"] || over["emax"]) {
      // A flat spectrum spanning the requested range, so the table covers exactly that range.
      YAML::Node f;
      f["type"] = "flat";
      f["emin"] = over["emin"] ? over["emin"].as<double>() : 0.1;
      f["emax"] = over["emax"] ? over["emax"].as<double>() : 20.0;
      config["beam"]["flux"] = f;
      config["xsec_table"]["emin"] = f["emin"];
    }
    if (over["nodes"]) config["xsec_table"]["nodes"] = over["nodes"];
    if (over["samples"]) config["xsec_table"]["samples"] = over["samples"];
    Generator gen(config);
    const auto& sel = gen.Selector();
    std::printf("# E[GeV]");
    for (size_t i = 0; i < sel.Size(); ++i)
      std::printf("  %s_%s[1e-38cm2/nucleus]", sel.Get(i).channel->Name().c_str(), sel.Get(i).nucleus.name.c_str());
    std::printf("  total[1e-38cm2/target]\n");
    const auto& e = sel.Get(0).table.energy;
    for (size_t j = 0; j < e.size(); ++j) {
      std::printf("%10.5f", e[j]);
      for (size_t i = 0; i < sel.Size(); ++i) std::printf("  %12.6g", sel.Get(i).table.sigma[j]);
      std::printf("  %12.6g\n", sel.TotalSigma(e[j]));
    }
    std::printf("# flux-averaged total: %.6g 1e-38 cm^2 per target unit\n", gen.Sampler().FluxAveragedSigma());
  } catch (const std::exception& ex) {
    std::fprintf(stderr, "nuflow-xsec: error: %s\n", ex.what());
    return 1;
  }
  return 0;
}
