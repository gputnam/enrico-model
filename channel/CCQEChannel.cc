#include "channel/CCQEChannel.h"

#include "core/PDG.h"

namespace nuflow {

// nu + n -> l- + p ; nubar + p -> l+ + n
bool CCQEChannel::Applies(int nuPdg, const Nucleus& nucleus) const {
  if (!pdg::IsNeutrino(nuPdg)) return false;
  return Multiplicity(nuPdg, nucleus) > 0;
}

HitNucleonRequest CCQEChannel::Request(int nuPdg) const {
  HitNucleonRequest r;
  r.nNucleons = 1;
  r.species[0] = nuPdg > 0 ? pdg::kNeutron : pdg::kProton;
  return r;
}

double CCQEChannel::Multiplicity(int nuPdg, const Nucleus& nucleus) const {
  return nuPdg > 0 ? nucleus.N() : nucleus.Z;
}

}  // namespace nuflow
