#include "channel/CCRESChannel.h"

#include "core/PDG.h"

namespace nuflow {

bool CCRESChannel::Applies(int nuPdg, const Nucleus& nucleus) const {
  return pdg::IsNeutrino(nuPdg) && nucleus.A >= 1;
}

HitNucleonRequest CCRESChannel::Request(int nuPdg) const {
  HitNucleonRequest r;
  r.nNucleons = 1;
  r.species[0] = pdg::kProton;
  return r;
}

HitNucleonRequest CCRESChannel::SampleRequest(int nuPdg, const Nucleus& nucleus, Random& rng) const {
  HitNucleonRequest r;
  r.nNucleons = 1;
  if (nucleus.Z == 0) r.species[0] = pdg::kNeutron;
  else if (nucleus.Z == nucleus.A) r.species[0] = pdg::kProton;
  else r.species[0] = rng.Uniform() * nucleus.A < nucleus.Z ? pdg::kProton : pdg::kNeutron;
  return r;
}

double CCRESChannel::Multiplicity(int nuPdg, const Nucleus& nucleus) const { return nucleus.A; }

}  // namespace nuflow
