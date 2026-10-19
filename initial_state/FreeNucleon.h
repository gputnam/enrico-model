#pragma once
// A free nucleon at rest (hydrogen, or a nucleon target): no binding, no Pauli blocking.

#include "initial_state/NuclearModel.h"

namespace nuflow {

class FreeNucleon : public NuclearModel {
 public:
  std::string Name() const override { return "FreeNucleon"; }
  void Sample(const Nucleus& nucleus, const HitNucleonRequest& request, InitialState& out,
              Random& rng) const override;
  bool PauliBlocked(const Nucleus&, int, const ThreeVector&) const override { return false; }
};

}  // namespace nuflow
