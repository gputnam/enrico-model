#pragma once
// Charged-current coherent pion production on a nucleus, nu A -> l- pi+ A(g.s.), nubar A -> l+ pi- A(g.s.)
// (Berger-Sehgal PCAC model, theory/coherent_pion.md). The nucleus acts as a whole: one interaction
// "site" per nucleus, no struck nucleon. Not defined on free nucleons (there the process is part of
// single-pion production, CCRES).

#include "channel/InteractionChannel.h"

namespace nuflow {

class CCCOHChannel : public InteractionChannel {
 public:
  using InteractionChannel::InteractionChannel;
  static constexpr int kProcessId = 100;  // NuHepMC: 100 = CC coherent

  std::string Name() const override { return "CCCOH"; }
  int ProcessId() const override { return kProcessId; }
  bool Applies(int nuPdg, const Nucleus& nucleus) const override;
  HitNucleonRequest Request(int nuPdg) const override;
  double Multiplicity(int nuPdg, const Nucleus& nucleus) const override { return 1.0; }
};

}  // namespace nuflow
