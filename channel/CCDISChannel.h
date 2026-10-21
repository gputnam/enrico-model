#pragma once
// Charged-current inclusive inelastic scattering above the resonance region (W > w_min), on either
// nucleon species: nu N -> l- X, nubar N -> l+ X, with X fragmented into a nucleon and pions.
// See theory/dis.md.

#include "channel/InteractionChannel.h"

namespace nuflow {

class CCDISChannel : public InteractionChannel {
 public:
  using InteractionChannel::InteractionChannel;
  static constexpr int kProcessId = 600;  // NuHepMC CC DIS range 600-649

  std::string Name() const override { return "CCDIS"; }
  int ProcessId() const override { return kProcessId; }
  bool Applies(int nuPdg, const Nucleus& nucleus) const override;
  HitNucleonRequest Request(int nuPdg) const override;
  // Proton with probability Z/A, neutron with probability N/A.
  HitNucleonRequest SampleRequest(int nuPdg, const Nucleus& nucleus, Random& rng) const override;
  double Multiplicity(int nuPdg, const Nucleus& nucleus) const override { return nucleus.A; }
};

}  // namespace nuflow
