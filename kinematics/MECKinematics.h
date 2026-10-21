#pragma once
// Charged-current 2p2h (meson-exchange current) kinematics on a nucleon pair from the stage-2 nuclear
// model:  nu + (np) -> l- + p p,  nu + (nn) -> l- + n p;  nubar + (pp) -> l+ + n p,  nubar + (np) -> l+ + n n.
// The fully differential cross section is evaluated event by event with the relativistic MEC of
// Ruiz Simo et al. (arXiv:1604.08423, see MECCurrent.h) in the Fermi gas: the lepton is sampled in
// (omega, |q|, phi), the outgoing pair in its rest frame, and the weight is the exact 2p2h integrand
// (direct and exchange, all spins) divided by the proposal density, so the nucleon-pair kinematics
// follow the model rather than an isotropic phase-space decay. See theory/mec.md.

#include "kinematics/KinematicsGenerator.h"
#include "kinematics/MECCurrent.h"

namespace nuflow {

class MECKinematics : public KinematicsGenerator {
 public:
  struct Options {
    double proposalQScale = 1.0;  // q0 of the |q| proposal ~ 1/(|q| + q0)^2 [GeV]
    double q3Max = 0;             // optional upper limit on |q| [GeV]; 0 = none
  };

  MECKinematics(MECCurrent current, Options opt) : current_(std::move(current)), opt_(opt) {}

  std::string Name() const override;
  double Generate(int nuPdg, double enu, const Nucleus& nucleus, const InitialState& initial,
                  const NuclearModel& nuclearModel, Event& event, Random& rng) const override;

  const MECCurrent& Current() const { return current_; }

 private:
  MECCurrent current_;
  Options opt_;
};

}  // namespace nuflow
