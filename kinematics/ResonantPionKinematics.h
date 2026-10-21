#pragma once
// Charged-current single-pion production in the Delta(1232) region on a nucleon from the stage-2
// nuclear model (impulse approximation, de Forest prescription for the off-shell nucleon):
//   nu p -> l- p pi+,  nu n -> l- n pi+,  nu n -> l- p pi0     (and nubar p/n -> l+ N pi-/pi0)
// with the HNV hadronic current (kinematics/PionProductionCurrent). The full three-body final state
// (lepton, nucleon, pion) is sampled, so the pion angular distribution follows from the matrix
// element. See theory/pion_production.md.

#include "kinematics/KinematicsGenerator.h"
#include "kinematics/PionProductionCurrent.h"

namespace nuflow {

class ResonantPionKinematics : public KinematicsGenerator {
 public:
  struct Options {
    double wMax = 1.6;             // upper limit of the generated pi N invariant mass [GeV]
    double proposalQ2Scale = 1.0;  // scale of the (1 + Q2/scale)^-2 proposal in the CM lepton angle [GeV^2]
    double bwFraction = 0.7;       // fraction of W proposals drawn from a Delta Breit-Wigner
  };

  ResonantPionKinematics(PionProductionCurrent current, Options opt) : current_(std::move(current)), opt_(opt) {}

  std::string Name() const override;
  double Generate(int nuPdg, double enu, const Nucleus& nucleus, const InitialState& initial,
                  const NuclearModel& nuclearModel, Event& event, Random& rng) const override;

  const PionProductionCurrent& Current() const { return current_; }

 private:
  PionProductionCurrent current_;
  Options opt_;
};

}  // namespace nuflow
