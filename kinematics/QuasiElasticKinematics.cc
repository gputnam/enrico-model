#include "kinematics/QuasiElasticKinematics.h"

#include <algorithm>
#include <cmath>

#include "core/Constants.h"
#include "core/PDG.h"
#include "kinematics/QECurrent.h"

namespace nuflow {

using namespace constants;

namespace {

// Orthonormal basis (e1, e2) perpendicular to unit vector n.
void Perpendicular(const ThreeVector& n, ThreeVector& e1, ThreeVector& e2) {
  const ThreeVector a = std::abs(n.x) < 0.9 ? ThreeVector{1, 0, 0} : ThreeVector{0, 1, 0};
  e1 = n.Cross(a).Unit();
  e2 = n.Cross(e1);
}

}  // namespace

double QuasiElasticKinematics::Generate(int nuPdg, double enu, const Nucleus& nucleus, const InitialState& initial,
                                        const NuclearModel& nuclearModel, Event& event, Random& rng) const {
  const HitNucleon& hit = initial.hit[0];
  const int helicity = nuPdg > 0 ? +1 : -1;
  const int leptonPdg = pdg::ChargedLeptonFor(nuPdg);
  const int finalPdg = nuPdg > 0 ? pdg::kProton : pdg::kNeutron;
  const double mLep = pdg::LeptonMass(leptonPdg);
  const double mIn = pdg::NucleonMass(hit.pdg);
  const double mOut = pdg::NucleonMass(finalPdg);

  // Two-body final state in the rest frame of k + p_bound (the bound nucleon is off shell).
  const FourVector k(enu, 0, 0, enu);
  const FourVector pBound = hit.Momentum();
  const FourVector Ptot = k + pBound;
  const double s = Ptot.M2();
  if (Ptot.e <= 0 || s <= (mLep + mOut) * (mLep + mOut)) return 0;
  const double sqrtS = std::sqrt(s);
  const ThreeVector beta = Ptot.BoostVector();

  const FourVector kCM = k.Boosted(-beta);
  const double kStar = kCM.e;
  const double eLepStar = (s + mLep * mLep - mOut * mOut) / (2 * sqrtS);
  const double pLepStar = std::sqrt(std::max(0.0, eLepStar * eLepStar - mLep * mLep));

  // Proposal for x = 1 - cos(theta*), theta* measured from the neutrino direction in the CM frame.
  // Q2 is linear in x (Q2 = Q2(x=0) + 2 k* p* x), so g(x) ~ (1 + b x)^-2 with b = 2 k* p* / scale.
  const double bCoef = 2 * kStar * pLepStar / opt_.proposalQ2Scale;
  double x, g;
  if (bCoef < 1e-6) {
    x = 2 * rng.Uniform();
    g = 0.5;
  } else {
    const double kappa = 2 * bCoef / (1 + 2 * bCoef);
    x = (1.0 / (1.0 - rng.Uniform() * kappa) - 1.0) / bCoef;
    const double d = 1 + bCoef * x;
    g = bCoef / (kappa * d * d);
  }
  const double cosT = std::clamp(1.0 - x, -1.0, 1.0);
  const double sinT = std::sqrt(std::max(0.0, 1 - cosT * cosT));
  const double phi = 2 * kPi * rng.Uniform();

  const ThreeVector n = kCM.Vect().Unit();
  ThreeVector e1, e2;
  Perpendicular(n, e1, e2);
  const ThreeVector dir = n * cosT + e1 * (sinT * std::cos(phi)) + e2 * (sinT * std::sin(phi));
  const FourVector lepton = FourVector(eLepStar, dir * pLepStar).Boosted(beta);
  const FourVector nucleonOut = Ptot - lepton;

  if (nuclearModel.PauliBlocked(nucleus, finalPdg, nucleonOut.Vect())) return 0;

  // Hadronic vertex: de Forest's prescription uses the on-shell initial nucleon and q~ = p' - p_on.
  const double eOn = hit.OnShellEnergy(mIn);
  const FourVector pOn(eOn, hit.p);
  const FourVector qTilde = nucleonOut - pOn;
  const double Q2Tilde = -qTilde.M2();
  if (qTilde.e <= 0 || Q2Tilde <= 0) return 0;

  const WeakFormFactors ff = ff_.Evaluate(Q2Tilde);
  const double M = 0.5 * (mIn + mOut);
  const double LH = ContractQE(k, lepton, pOn, nucleonOut, ff, M, helicity);
  if (LH <= 0) return 0;

  // dsigma/dOmega* = <|M|^2> / (4 E_nu E_N) * |p*| / (16 pi^2 sqrt(s)),  <|M|^2> = G_F^2 Vud^2 LH / 4
  const double dSigmaDOmega =
      kGF * kGF * kVud * kVud * 0.25 * LH / (4 * enu * eOn) * pLepStar / (16 * kPi * kPi * sqrtS);
  const double weight = dSigmaDOmega * (2 * kPi) / g * kGeV2To1e38Cm2;

  // Fill the primary vertex.
  const FourVector q = k - lepton;
  event.nuPdg = nuPdg;
  event.targetPdg = nucleus.pdg;
  event.enu = enu;
  event.Q2 = -q.M2();
  event.q0 = q.e;
  event.q3 = q.P();
  event.W = std::sqrt(std::max(0.0, (pOn + qTilde).M2()));
  event.initial = initial;
  auto& parts = event.particles;
  parts.clear();
  event.Add(nuPdg, kStatusBeam, k);
  event.Add(nucleus.pdg, kStatusTarget, FourVector(nucleus.mass, 0, 0, 0));
  const int iHit = event.Add(hit.pdg, kStatusStruckNucleon, pBound, 1);
  event.Add(leptonPdg, kStatusFinal, lepton, 0);
  event.Add(finalPdg, kStatusPreFSI, nucleonOut, iHit);
  if (initial.remnantPdg != 0) event.Add(initial.remnantPdg, kStatusRemnant, initial.remnant, 1);
  return weight;
}

}  // namespace nuflow
