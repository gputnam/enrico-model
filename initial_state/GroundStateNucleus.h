#pragma once
// Stage 2 for coherent processes: the whole nucleus, at rest and in its ground state, takes part in
// the interaction; no nucleon is singled out (nHit = 0). The remnant is the full target nucleus.

#include "initial_state/NuclearModel.h"

namespace nuflow {

class GroundStateNucleus : public NuclearModel {
 public:
  std::string Name() const override { return "GroundStateNucleus"; }
  void Sample(const Nucleus& nucleus, const HitNucleonRequest& request, InitialState& out,
              Random& rng) const override;
  bool PauliBlocked(const Nucleus&, int, const ThreeVector&) const override { return false; }
};

}  // namespace nuflow
