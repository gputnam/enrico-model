#pragma once
// Charged-current inclusive inelastic scattering (shallow- and deep-inelastic, W > w_min) on a nucleon
// from the stage-2 nuclear model, from the Bodek-Yang structure functions (kinematics/BodekYang),
// with the de Forest prescription for the off-shell nucleon (as in CCQE and CCRES). The lepton is
// sampled in the CM frame of k + p_bound against a hadronic system of mass W; after acceptance the
// hadronic system is fragmented into a nucleon and pions (kinematics/DISHadronization).
// See theory/dis.md.

#include "kinematics/BodekYang.h"
#include "kinematics/DISHadronization.h"
#include "kinematics/KinematicsGenerator.h"

namespace nuflow {

class DISKinematics : public KinematicsGenerator {
 public:
  struct Options {
    double wMin = 1.6;             // lower limit of the hadronic invariant mass [GeV] (= CCRES w_max)
    double proposalQ2Scale = 1.0;  // Q2 scale of the peaked part of the lepton-angle proposal [GeV^2]
    double flatAngleFraction = 0.5;  // fraction of lepton-angle proposals flat in cos(theta*)
  };

  DISKinematics(BodekYang sf, DISHadronization had, Options opt)
      : sf_(std::move(sf)), had_(std::move(had)), opt_(opt) {}

  std::string Name() const override;
  double Generate(int nuPdg, double enu, const Nucleus& nucleus, const InitialState& initial,
                  const NuclearModel& nuclearModel, Event& event, Random& rng) const override;
  // Fragments the hadronic system of an accepted event (appends status-22 hadrons).
  void Finalize(const Nucleus& nucleus, Event& event, Random& rng) const override;

  // L_{mu nu} W^{mu nu} for the general hadronic tensor with structure functions F1..F5 of a
  // nucleon with four-momentum P (mass M) and momentum transfer q; k, kp: neutrino and charged lepton.
  // helicity = +1 (nu), -1 (nubar). The rest-frame cross section is
  //   d2sigma/dx dy = G_F^2 M nu / (16 pi E) * Contract(...)  (checked in tests/test_dis.cc).
  static double Contract(const FourVector& k, const FourVector& kp, const FourVector& P, double M,
                         const FourVector& q, const StructureFunctions& F, int helicity);

  const BodekYang& StructureFunctionModel() const { return sf_; }
  const DISHadronization& Hadronization() const { return had_; }

 private:
  BodekYang sf_;
  DISHadronization had_;
  Options opt_;
};

}  // namespace nuflow
