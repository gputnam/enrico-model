#include "channel/CCCOHChannel.h"

#include "core/PDG.h"

namespace nuflow {

bool CCCOHChannel::Applies(int nuPdg, const Nucleus& nucleus) const {
  return pdg::IsNeutrino(nuPdg) && nucleus.A >= 4;
}

HitNucleonRequest CCCOHChannel::Request(int) const {
  HitNucleonRequest r;
  r.nNucleons = 0;
  return r;
}

}  // namespace nuflow
