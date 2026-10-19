#pragma once
// Orchestration of the four generation stages:
//   E_nu ~ phi(E) sigma_tot(E)  ->  stage 1: (channel, nucleus) ~ n_A sigma_{c,A}(E)
//   -> stage 2: initial state from the nuclear model  -> stage 3: kinematics, accept w.p. w / w_max(E)
//   (stages 2-3 are repeated at fixed E and channel until accepted)  -> stage 4: FSI.
// Events are produced in fixed-size batches, each with its own random stream, and written in batch
// order, so output is identical for a given seed regardless of the number of threads.

#include <memory>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "channel/ChannelSelector.h"
#include "core/Nucleus.h"
#include "flux/FluxSource.h"
#include "fsi/FSIModel.h"
#include "orchestration/EnergySampler.h"
#include "orchestration/OutputWriter.h"

namespace nuflow {

class Generator {
 public:
  // `config` is a full run configuration (see config/*.yaml).
  explicit Generator(const YAML::Node& config);
  ~Generator();

  // Generate `nEvents` events and hand them to the writers.
  RunSummary Run(long long nEvents, const std::vector<OutputWriter*>& writers);

  const ChannelSelector& Selector() const { return *selector_; }
  const EnergySampler& Sampler() const { return *sampler_; }
  RunInfo Info() const;

  // Generate a single event with the given random stream (exposed for tests).
  // Returns the number of stage-2/3 trials used. Envelope violations are counted in nOverweight and
  // the largest w / w_envelope seen in any trial is tracked in maxRatio (> 1 means a violation).
  long GenerateOne(Random& rng, Event& event, long long* nOverweight = nullptr, double* maxRatio = nullptr) const;

 private:
  YAML::Node config_;
  int nuPdg_ = 0;
  uint64_t seed_ = 1;
  int threads_ = 1;
  int batchSize_ = 5000;
  double safetyFactor_ = 1.2;
  Target target_;
  std::unique_ptr<FluxSource> flux_;
  std::vector<std::shared_ptr<InteractionChannel>> channels_;
  std::unique_ptr<ChannelSelector> selector_;
  std::unique_ptr<EnergySampler> sampler_;
  std::unique_ptr<FSIModel> fsi_;
};

}  // namespace nuflow
