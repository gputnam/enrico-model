#pragma once
// Stage 3: kinematics generation.
//
// Given the neutrino, the nucleus and the stage-2 initial state, a KinematicsGenerator samples the
// final-state kinematics of its channel from a proposal density and returns the corresponding
// event weight, i.e. the differential cross section divided by the proposal density. The weight is
// in units of 1e-38 cm^2 per struck nucleon (or nucleon pair); its expectation over the stage-2
// distribution is the total cross section per struck nucleon. A weight of zero means the trial is
// forbidden (below threshold, Pauli blocked, ...).

#include <string>

#include "core/Event.h"
#include "core/Nucleus.h"
#include "core/Random.h"
#include "initial_state/NuclearModel.h"

namespace nuflow {

class KinematicsGenerator {
 public:
  virtual ~KinematicsGenerator() = default;
  virtual std::string Name() const = 0;

  // Fills `event.particles` with the primary-vertex particles and the event-level kinematic summary.
  virtual double Generate(int nuPdg, double enu, const Nucleus& nucleus, const InitialState& initial,
                          const NuclearModel& nuclearModel, Event& event, Random& rng) const = 0;

  // Called once for each *accepted* event, before stage 4: completes the primary-vertex final state
  // where that is independent of the weight (e.g. DIS hadronization). The default does nothing and
  // consumes no random numbers.
  virtual void Finalize(const Nucleus& nucleus, Event& event, Random& rng) const {}
};

}  // namespace nuflow
