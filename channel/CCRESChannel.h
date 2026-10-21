#pragma once
// Charged-current single-pion production in the Delta(1232) region (resonant + chiral background),
// on either nucleon species:  nu p -> l- p pi+,  nu n -> l- n pi+ / l- p pi0  (nubar: l+ N pi-/pi0).
// See theory/pion_production.md.

#include "channel/InteractionChannel.h"

namespace nuflow {

class CCRESChannel : public InteractionChannel {
 public:
  using InteractionChannel::InteractionChannel;
  static constexpr int kProcessId = 400;  // NuHepMC CC resonant pion production range 400-449

  std::string Name() const override { return "CCRES"; }
  int ProcessId() const override { return kProcessId; }
  bool Applies(int nuPdg, const Nucleus& nucleus) const override;
  HitNucleonRequest Request(int nuPdg) const override;
  // Proton with probability Z/A, neutron with probability N/A.
  HitNucleonRequest SampleRequest(int nuPdg, const Nucleus& nucleus, Random& rng) const override;
  // All A nucleons are candidates; the per-trial weight averages over the species draw.
  double Multiplicity(int nuPdg, const Nucleus& nucleus) const override;
};

}  // namespace nuflow
