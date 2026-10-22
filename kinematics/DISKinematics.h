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
    double wMax = 0;               // upper limit of W [GeV]; 0 = kinematic limit
    // Exclusive-channel decomposition of the inelastic continuum (theory/dis.md section 3). S(W) = P_1pi(W)
    // x (Bodek-Yang strength) is the single-pion (N pi) part implied by the hadronization model; it is split
    // into isospin components using the pure-I=3/2 partner target (nu p, nubar n) at the same kinematics,
    // S_3 = S_partner/3 on the mixed target (nu n, nubar p) and S_1 = S - S_3, and each is scaled:
    //   sigma_1pi = R_3/2 S_3 + R_1/2 S_1  (pure target: R_3/2 S).
    // The Clebsch-Gordan coefficients fix the charge states (I = 3/2: n pi+ : p pi0 = 1 : 2; I = 1/2: 2 : 1).
    // The rest of the Bodek-Yang strength goes to two or more pions (inclusive strength kept).
    //  - isospinDecomposition = false: legacy, sigma_1pi = R_3/2 S on every target, charges from the
    //    hadronization model.
    double scaleI32 = 1.0;
    double scaleI12 = 1.0;
    bool isospinDecomposition = false;
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

  // Single-pion probability among the final states of the Bodek-Yang strength at this kinematic point,
  // and the probability of the charged-pion state (nu n -> n pi+, nubar p -> p pi-) among them.
  // lwSelf, lwPartner: L.W of the struck nucleon and of its isospin partner (same kinematics).
  struct SinglePionSplit {
    double p1 = 0;
    double chargedShare = 1;
  };
  SinglePionSplit SinglePion(int nuPdg, int hitPdg, double W, double lwSelf, double lwPartner) const;
  // True if SinglePion needs the partner's L.W for this neutrino and nucleon.
  bool NeedsPartner(int nuPdg, int hitPdg) const;

  const BodekYang& StructureFunctionModel() const { return sf_; }
  const DISHadronization& Hadronization() const { return had_; }

 private:
  // L.W for the nucleon `proton` at the kinematics of a struck nucleon with on-shell momentum pOn.
  double ContractFor(bool neutrino, bool proton, const FourVector& k, const FourVector& lepton, const FourVector& pOn,
                     double mIn, double W) const;
  BodekYang sf_;
  DISHadronization had_;
  Options opt_;
};

}  // namespace nuflow
