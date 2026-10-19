#pragma once
// Charged-current quasi-elastic scattering on a single nucleon.

#include "channel/InteractionChannel.h"

namespace nuflow {

class CCQEChannel : public InteractionChannel {
 public:
  using InteractionChannel::InteractionChannel;
  static constexpr int kProcessId = 200;  // NuHepMC CCQE range 200-249

  std::string Name() const override { return "CCQE"; }
  int ProcessId() const override { return kProcessId; }
  bool Applies(int nuPdg, const Nucleus& nucleus) const override;
  HitNucleonRequest Request(int nuPdg) const override;
  double Multiplicity(int nuPdg, const Nucleus& nucleus) const override;
};

}  // namespace nuflow
