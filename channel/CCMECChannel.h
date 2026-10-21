#pragma once
// Charged-current two-particle-two-hole (2p2h) excitation by meson-exchange currents, on a nucleon pair:
//   nu:    (n p) -> l- p p,   (n n) -> l- n p
//   nubar: (p p) -> l+ n p,   (n p) -> l+ n n
// Relativistic MEC of Ruiz Simo et al. (arXiv:1604.08423) in the Fermi gas. See theory/mec.md.

#include "channel/InteractionChannel.h"

namespace nuflow {

class CCMECChannel : public InteractionChannel {
 public:
  using InteractionChannel::InteractionChannel;
  static constexpr int kProcessId = 300;  // NuHepMC CC MEC/2p2h range 300-349

  std::string Name() const override { return "CCMEC"; }
  int ProcessId() const override { return kProcessId; }
  // Nuclei with A >= 3 and a Fermi sea (the free nucleon and the deuteron have no 2p2h here).
  bool Applies(int nuPdg, const Nucleus& nucleus) const override;
  HitNucleonRequest Request(int nuPdg) const override;
  // Draws the pair type with probability 1/2 each: nu (n,p) or (n,n); nubar (n,p) or (p,p).
  HitNucleonRequest SampleRequest(int nuPdg, const Nucleus& nucleus, Random& rng) const override;
  // 2N (nu) or 2Z (nubar): with the 1/2 pair-type draw, Multiplicity x <w> = N (w_np + w_nn) (theory/mec.md).
  double Multiplicity(int nuPdg, const Nucleus& nucleus) const override;
};

}  // namespace nuflow
