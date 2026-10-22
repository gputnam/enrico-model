#include "channel/CCNStarChannel.h"

#include "core/PDG.h"

namespace nuflow {

bool CCNStarChannel::Applies(int nuPdg, const Nucleus& nucleus) const {
  return pdg::IsNeutrino(nuPdg) && Multiplicity(nuPdg, nucleus) > 0;
}

HitNucleonRequest CCNStarChannel::Request(int nuPdg) const {
  HitNucleonRequest r;
  r.nNucleons = 1;
  r.species[0] = nuPdg > 0 ? pdg::kNeutron : pdg::kProton;
  return r;
}

double CCNStarChannel::Multiplicity(int nuPdg, const Nucleus& nucleus) const {
  return nuPdg > 0 ? nucleus.A - nucleus.Z : nucleus.Z;
}

}  // namespace nuflow
