#pragma once
// Charged-current quasi-elastic kinematics:  nu + n -> l- + p,  nubar + p -> l+ + n,
// on a nucleon from the stage-2 nuclear model (impulse approximation, de Forest prescription).
// See theory/quasielastic.md.

#include "kinematics/FormFactors.h"
#include "kinematics/KinematicsGenerator.h"

namespace nuflow {

class QuasiElasticKinematics : public KinematicsGenerator {
 public:
  struct Options {
    // Scale [GeV^2] of the proposal density in 1 - cos(theta*) (CM lepton angle). The proposal is
    // ~ (1 + Q2/scale)^-2, flat for small Q2 and broader than the dipole fall-off of dsigma/dQ2.
    double proposalQ2Scale = 0.6;
  };

  QuasiElasticKinematics(NucleonFormFactors ff, Options opt) : ff_(std::move(ff)), opt_(opt) {}
  explicit QuasiElasticKinematics(NucleonFormFactors ff) : QuasiElasticKinematics(std::move(ff), Options{}) {}

  std::string Name() const override { return "CCQE[" + ff_.Name() + "]"; }
  double Generate(int nuPdg, double enu, const Nucleus& nucleus, const InitialState& initial,
                  const NuclearModel& nuclearModel, Event& event, Random& rng) const override;

  const NucleonFormFactors& FormFactors() const { return ff_; }

 private:
  NucleonFormFactors ff_;
  Options opt_;
};

}  // namespace nuflow
