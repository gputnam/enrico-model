#pragma once
// Stage 1 sampling: choose (channel, nucleus) at a given neutrino energy with probability
// proportional to  n_A * sigma_{c,A}(E),  from tabulated total cross sections.

#include <string>
#include <vector>

#include "channel/InteractionChannel.h"
#include "core/Nucleus.h"
#include "core/Random.h"

namespace nuflow {

// sigma(E) and maximum trial weight w_max(E) on an energy grid, for one (channel, nucleus) pair.
struct XSecTable {
  std::vector<double> energy;    // [GeV], increasing
  std::vector<double> sigma;     // per nucleus [1e-38 cm^2]
  std::vector<double> sigmaErr;  // MC integration uncertainty [1e-38 cm^2]
  std::vector<double> wMax;      // maximum per-nucleon trial weight seen at the node [1e-38 cm^2]

  double Sigma(double e) const;
  // Envelope for accept/reject between nodes: larger of the bracketing node maxima.
  double WMax(double e) const;

 private:
  size_t Bracket(double e) const;
};

// A (channel, nucleus) pair that can occur for the configured neutrino and target.
struct Process {
  const InteractionChannel* channel = nullptr;
  Nucleus nucleus;
  double count = 1;         // abundance of this nucleus in the target unit
  double multiplicity = 1;  // struck-nucleon multiplicity (e.g. N for nu CCQE)
  std::string physicsKey;   // full physics configuration (keys the cross-section cache)
  XSecTable table;
};

class ChannelSelector {
 public:
  explicit ChannelSelector(std::vector<Process> processes) : processes_(std::move(processes)) {}

  size_t Size() const { return processes_.size(); }
  const Process& Get(size_t i) const { return processes_[i]; }
  std::vector<Process>& Processes() { return processes_; }

  // Total cross section per target unit (e.g. per CH) [1e-38 cm^2].
  double TotalSigma(double e) const;
  // Index of the selected process (requires TotalSigma(e) > 0).
  size_t Select(double e, Random& rng) const;

 private:
  std::vector<Process> processes_;
};

}  // namespace nuflow
