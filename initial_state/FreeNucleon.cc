#include "initial_state/FreeNucleon.h"

#include <stdexcept>

#include "core/PDG.h"

namespace nuflow {

void FreeNucleon::Sample(const Nucleus& nucleus, const HitNucleonRequest& request, InitialState& out,
                         Random&) const {
  if (request.nNucleons != 1 || nucleus.A != 1)
    throw std::logic_error("FreeNucleon: only single-nucleon interactions on A=1 targets");
  out.nHit = 1;
  auto& h = out.hit[0];
  h.pdg = request.species[0];
  h.p = {0, 0, 0};
  h.energy = pdg::NucleonMass(h.pdg);
  h.removalEnergy = 0;
  out.remnantPdg = 0;
  out.remnant = {};
}

}  // namespace nuflow
