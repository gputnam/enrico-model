#include "initial_state/GroundStateNucleus.h"

namespace nuflow {

void GroundStateNucleus::Sample(const Nucleus& nucleus, const HitNucleonRequest&, InitialState& out, Random&) const {
  out = InitialState{};
  out.nHit = 0;
  out.remnantPdg = nucleus.pdg;
  out.remnant = FourVector(nucleus.mass, 0, 0, 0);
}

}  // namespace nuflow
