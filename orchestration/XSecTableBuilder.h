#pragma once
// Builds the per-(channel, nucleus) total cross-section tables by Monte Carlo integration of the
// SAME stage-2 x stage-3 trial weight used for event generation, so totals and differential
// distributions are consistent by construction. Tables are cached on disk keyed by a hash of the
// physics configuration.

#include <string>
#include <vector>

#include "channel/ChannelSelector.h"

namespace nuflow {

struct XSecTableOptions {
  std::vector<double> energies;  // grid [GeV]
  long samplesPerNode = 200000;
  uint64_t seed = 1;
  int threads = 1;
  std::string cacheDir;  // empty = no caching
  bool verbose = true;
};

// Log-spaced energy grid from emin to emax (n nodes), or a single node for a mono-energetic beam.
std::vector<double> MakeEnergyGrid(double emin, double emax, int n);

void BuildXSecTables(int nuPdg, std::vector<Process>& processes, const XSecTableOptions& opt);

}  // namespace nuflow
