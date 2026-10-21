#include "channel/CCMECChannel.h"

#include "core/PDG.h"

namespace nuflow {

bool CCMECChannel::Applies(int nuPdg, const Nucleus& nucleus) const {
  return pdg::IsNeutrino(nuPdg) && nucleus.A >= 3 && nucleus.kFermiProton > 0 && nucleus.kFermiNeutron > 0 &&
         Multiplicity(nuPdg, nucleus) > 0;
}

HitNucleonRequest CCMECChannel::Request(int nuPdg) const {
  HitNucleonRequest r;
  r.nNucleons = 2;
  r.species[0] = pdg::kNeutron;
  r.species[1] = pdg::kProton;
  return r;
}

HitNucleonRequest CCMECChannel::SampleRequest(int nuPdg, const Nucleus& nucleus, Random& rng) const {
  HitNucleonRequest r = Request(nuPdg);
  if (rng.Uniform() < 0.5) {
    const int same = nuPdg > 0 ? pdg::kNeutron : pdg::kProton;
    r.species[0] = r.species[1] = same;
  }
  return r;
}

double CCMECChannel::Multiplicity(int nuPdg, const Nucleus& nucleus) const {
  return nuPdg > 0 ? 2.0 * nucleus.N() : 2.0 * nucleus.Z;
}

}  // namespace nuflow
