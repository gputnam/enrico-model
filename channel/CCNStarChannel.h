#pragma once
// Charged-current excitation of the isospin-1/2 resonances N(1440), N(1520), N(1535) (second resonance
// region): nu n -> l- R+, nubar p -> l+ R0, followed by R -> N pi, N eta, N pi pi. The W+ (W-) cannot
// excite an I = 1/2 state on a proton (neutron). See theory/nstar.md.

#include "channel/InteractionChannel.h"

namespace nuflow {

class CCNStarChannel : public InteractionChannel {
 public:
  using InteractionChannel::InteractionChannel;
  static constexpr int kProcessId = 401;  // NuHepMC CC resonant range 400-449 (400 = CCRES Delta region)

  std::string Name() const override { return "CCNSTAR"; }
  int ProcessId() const override { return kProcessId; }
  bool Applies(int nuPdg, const Nucleus& nucleus) const override;
  // Neutron for neutrinos, proton for antineutrinos.
  HitNucleonRequest Request(int nuPdg) const override;
  // Number of neutrons (nu) or protons (nubar).
  double Multiplicity(int nuPdg, const Nucleus& nucleus) const override;
};

}  // namespace nuflow
