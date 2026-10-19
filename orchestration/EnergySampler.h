#pragma once
// Samples the neutrino energy of interacting neutrinos, distributed as  phi(E) * sigma_tot(E),
// where sigma_tot is the target-unit total cross section summed over all processes.

#include <vector>

#include "channel/ChannelSelector.h"
#include "core/Random.h"
#include "flux/FluxSource.h"

namespace nuflow {

class EnergySampler {
 public:
  EnergySampler(const FluxSource& flux, const ChannelSelector& selector);

  double Sample(Random& rng) const;

  // Flux-averaged total cross section per target unit [1e-38 cm^2].
  double FluxAveragedSigma() const { return fluxAveragedSigma_; }

 private:
  const ChannelSelector& selector_;
  bool mono_ = false;
  double monoEnergy_ = 0;
  std::vector<double> edges_;
  std::vector<double> cdf_;       // cumulative of  phi_i * int_bin sigma(E) dE
  std::vector<double> sigmaMax_;  // max of sigma in each bin (piecewise-linear => at nodes/edges)
  double fluxAveragedSigma_ = 0;
};

}  // namespace nuflow
