#include "kinematics/DISKinematics.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "core/Constants.h"
#include "core/PDG.h"

namespace nuflow {

using namespace constants;

namespace {

void Perpendicular(const ThreeVector& n, ThreeVector& e1, ThreeVector& e2) {
  const ThreeVector a = std::abs(n.x) < 0.9 ? ThreeVector{1, 0, 0} : ThreeVector{0, 1, 0};
  e1 = n.Cross(a).Unit();
  e2 = n.Cross(e1);
}

double Kallen(double a, double b, double c) { return a * a + b * b + c * c - 2 * (a * b + a * c + b * c); }

}  // namespace

std::string DISKinematics::Name() const {
  std::ostringstream s;
  s.precision(6);
  s << "CCDIS[" << sf_.Name() << ",Wmin=" << opt_.wMin;
  if (opt_.wMax > 0) s << ",Wmax=" << opt_.wMax;
  if (opt_.scaleI32 != 1.0 || opt_.scaleI12 != 1.0 || opt_.isospinDecomposition)
    s << ",R1pi=" << opt_.scaleI32 << "/" << opt_.scaleI12 << (opt_.isospinDecomposition ? "iso" : "");
  s << "]";
  return s.str();
}

bool DISKinematics::NeedsPartner(int nuPdg, int hitPdg) const {
  const bool mixed = (nuPdg > 0) == (hitPdg == pdg::kNeutron);  // nu n or nubar p
  return opt_.isospinDecomposition && mixed;
}

DISKinematics::SinglePionSplit DISKinematics::SinglePion(int nuPdg, int hitPdg, double W, double lwSelf,
                                                         double lwPartner) const {
  SinglePionSplit r;
  if (!(lwSelf > 0)) return r;
  const double sSelf = had_.SinglePionFraction(nuPdg, hitPdg, W) * lwSelf;
  if (!opt_.isospinDecomposition) {
    r.p1 = std::min(1.0, opt_.scaleI32 * sSelf / lwSelf);
    return r;
  }
  if (!NeedsPartner(nuPdg, hitPdg)) {  // nu p, nubar n: pure I = 3/2
    r.p1 = std::min(1.0, opt_.scaleI32 * sSelf / lwSelf);
    return r;
  }
  // Mixed target. nu p -> p pi+ (amplitude A3) is pure I = 3/2; on the neutron the I = 3/2 part gives
  // |A3|^2/9 in n pi+ and 2|A3|^2/9 in p pi0 (1/3 of the proton in total); the I = 1/2 remainder goes
  // 2/3 to n pi+ and 1/3 to p pi0 (Clebsch-Gordan; I = 1/2 - 3/2 interference neglected). Same for nubar
  // with p <-> n and pi+ -> pi-.
  const int partnerPdg = hitPdg == pdg::kProton ? pdg::kNeutron : pdg::kProton;
  const double sPartner = had_.SinglePionFraction(nuPdg, partnerPdg, W) * std::max(0.0, lwPartner);
  const double s3 = std::min(sSelf, sPartner / 3), s1 = sSelf - s3;
  const double t3 = opt_.scaleI32 * s3, t1 = opt_.scaleI12 * s1, tot = t3 + t1;
  r.p1 = std::min(1.0, tot / lwSelf);
  r.chargedShare = tot > 0 ? (t3 / 3 + 2 * t1 / 3) / tot : 1.0 / 3;
  return r;
}

double DISKinematics::ContractFor(bool neutrino, bool proton, const FourVector& k, const FourVector& lepton,
                                  const FourVector& pOn, double mIn, double W) const {
  const FourVector q = k - lepton;
  const double Q2 = -q.M2(), w2m2 = W * W - mIn * mIn + Q2;
  if (!(Q2 > 0) || !(w2m2 > 0)) return 0;
  const double xEff = Q2 / w2m2, nuEff = w2m2 / (2 * mIn);
  return Contract(k, lepton, pOn, mIn, q, sf_.Compute(neutrino, proton, xEff, Q2, nuEff, mIn), neutrino ? +1 : -1);
}

double DISKinematics::Contract(const FourVector& k, const FourVector& kp, const FourVector& P, double M,
                               const FourVector& q, const StructureFunctions& F, int helicity) {
  // W^{mu nu} = -g W1 + P P/M^2 W2 - i eps^{mu nu a b} P_a q_b/(2M^2) W3 + q q/M^2 W4 + (P q + q P)/(2M^2) W5
  // with W1 = F1/M, W2 = F2/nu, W3 = F3/nu, W4 = F4/(x nu), W5 = 2 F5/nu (nu = P.q/M, x = -q^2/2P.q):
  // the normalisation of F4, F5 in which the rest-frame d2sigma/dxdy has the form of Kretzer & Reno
  // (PRD 66, 113007) and the LO (Albright-Jarlskog) relations read F4 = 0, 2x F5 = F2 (W5 = W2/x).
  // L^{mu nu} = 8 [k k' + k' k - g (k.k') -/+ i eps k k'] (the lepton mass drops out of the trace).
  const double nu = P.Dot(q) / M;
  if (!(nu > 0)) return 0;
  const double xB = -q.M2() / (2 * P.Dot(q));
  const double W1 = F.F1 / M, W2 = F.F2 / nu, W3 = F.F3 / nu, W4 = xB > 0 ? F.F4 / (xB * nu) : 0, W5 = 2 * F.F5 / nu;
  const double kkp = k.Dot(kp), kP = k.Dot(P), kpP = kp.Dot(P), kq = k.Dot(q), kpq = kp.Dot(q), Pq = P.Dot(q),
               q2 = q.M2(), M2 = M * M;
  const double T1 = 16 * kkp;
  const double T2 = 8 * (2 * kP * kpP - M2 * kkp) / M2;
  const double T3 = 8 * helicity * (kP * kpq - kq * kpP) / M2;
  const double T4 = 8 * (2 * kq * kpq - q2 * kkp) / M2;
  const double T5 = 8 * (kP * kpq + kq * kpP - Pq * kkp) / M2;
  return W1 * T1 + W2 * T2 + W3 * T3 + W4 * T4 + W5 * T5;
}

double DISKinematics::Generate(int nuPdg, double enu, const Nucleus& nucleus, const InitialState& initial,
                               const NuclearModel& nuclearModel, Event& event, Random& rng) const {
  const HitNucleon& hit = initial.hit[0];
  const bool nu = nuPdg > 0;
  const int leptonPdg = pdg::ChargedLeptonFor(nuPdg);
  const double mLep = pdg::LeptonMass(leptonPdg);
  const double mIn = pdg::NucleonMass(hit.pdg);

  const FourVector k(enu, 0, 0, enu);
  const FourVector pBound = hit.Momentum();
  const FourVector Ptot = k + pBound;
  const double s = Ptot.M2();
  if (Ptot.e <= 0 || s <= 0) return 0;
  const double sqrtS = std::sqrt(s);
  const double wLo = std::max(opt_.wMin, kMassProton + kMassPionNeutral);
  const double wHi = opt_.wMax > 0 ? std::min(opt_.wMax, sqrtS - mLep) : sqrtS - mLep;
  if (wHi <= wLo) return 0;

  // --- W proposal: half flat in W, half flat in W^2.
  double W;
  if (rng.Uniform() < 0.5) W = wLo + rng.Uniform() * (wHi - wLo);
  else W = std::sqrt(wLo * wLo + rng.Uniform() * (wHi * wHi - wLo * wLo));
  const double gW = 0.5 / (wHi - wLo) + 0.5 * 2 * W / (wHi * wHi - wLo * wLo);

  // --- Lepton in the CM frame of k + p_bound.
  const ThreeVector beta = Ptot.BoostVector();
  const FourVector kCM = k.Boosted(-beta);
  const double kStar = kCM.e;
  const double lamL = Kallen(s, mLep * mLep, W * W);
  if (lamL <= 0) return 0;
  const double pLepStar = std::sqrt(lamL) / (2 * sqrtS);
  const double eLepStar = std::sqrt(pLepStar * pLepStar + mLep * mLep);

  // x = 1 - cos(theta*) in [0, 2]: mixture of flat and (1 + b x)^-2 (Q2 ~ 2 k* p* x).
  const double bCoef = 2 * kStar * pLepStar / opt_.proposalQ2Scale;
  const double f = opt_.flatAngleFraction;
  double x;
  const bool peaked = bCoef > 1e-6;
  const double kappa = peaked ? 2 * bCoef / (1 + 2 * bCoef) : 0;
  if (!peaked || rng.Uniform() < f) {
    x = 2 * rng.Uniform();
  } else {
    x = (1.0 / (1.0 - rng.Uniform() * kappa) - 1.0) / bCoef;
  }
  double gx = 0.5;
  if (peaked) {
    const double d = 1 + bCoef * x;
    gx = f * 0.5 + (1 - f) * bCoef / (kappa * d * d);
  }
  const double cosT = std::clamp(1.0 - x, -1.0, 1.0);
  const double sinT = std::sqrt(std::max(0.0, 1 - cosT * cosT));
  const double phi = 2 * kPi * rng.Uniform();
  const ThreeVector n = kCM.Vect().Unit();
  ThreeVector e1, e2;
  Perpendicular(n, e1, e2);
  const ThreeVector dir = n * cosT + e1 * (sinT * std::cos(phi)) + e2 * (sinT * std::sin(phi));
  const FourVector lepton = FourVector(eLepStar, dir * pLepStar).Boosted(beta);
  const FourVector q = k - lepton;
  const FourVector X = Ptot - lepton;  // hadronic system, mass W

  // --- Bound nucleon. The structure functions are those of a free nucleon at the same Q2 (lepton
  // side) and the same, physical, hadronic mass W (energy conservation with the off-shell nucleon):
  //   x_eff = Q2 / (W^2 - M^2 + Q2),  nu_eff = (W^2 - M^2 + Q2) / 2M.
  // The binding energy thus lowers W at fixed lepton kinematics. (A de Forest q~ = q - (eps, 0) would
  // shift Q2 by 2 eps q0, i.e. by ~1 GeV^2 at q0 ~ 20 GeV, which is unphysical for DIS.) The tensor
  // is built from the on-shell four-momentum p_on and q, so current conservation removes the q^mu
  // terms up to O(m_l^2), and Fermi motion enters through p_on and the flux factor.
  const double eOn = hit.OnShellEnergy(mIn);
  const FourVector pOn(eOn, hit.p);
  const double Q2 = -q.M2();
  if (!(Q2 > 0)) return 0;
  const double w2m2 = W * W - mIn * mIn + Q2;
  if (!(w2m2 > 0)) return 0;
  const double xEff = Q2 / w2m2, nuEff = w2m2 / (2 * mIn);
  const StructureFunctions F = sf_.Compute(nu, hit.pdg == pdg::kProton, xEff, Q2, nuEff, mIn);
  const double LW = Contract(k, lepton, pOn, mIn, q, F, nu ? +1 : -1);
  if (!(LW > 0)) return 0;

  // dsigma = G_F^2 M / (32 pi^2 E_nu E_N) L.W d^3k'/E' (propagator), d^3k'/E' = p*/(2 sqrt s) dW^2 dOmega*.
  const double prop = 1 / ((1 + Q2 / (kMassW * kMassW)) * (1 + Q2 / (kMassW * kMassW)));
  const double dsig = kGF * kGF * mIn / (32 * kPi * kPi * enu * eOn) * LW * prop * pLepStar / (2 * sqrtS) * 2 * W;
  const double weight = dsig * (2 * kPi) / (gx * gW) * kGeV2To1e38Cm2;

  // --- Primary vertex (hadrons are added in Finalize for accepted events).
  event.nuPdg = nuPdg;
  event.targetPdg = nucleus.pdg;
  event.enu = enu;
  event.Q2 = Q2;
  event.q0 = q.e;
  event.q3 = q.P();
  event.W = W;
  event.initial = initial;
  auto& parts = event.particles;
  parts.clear();
  event.Add(nuPdg, kStatusBeam, k);
  event.Add(nucleus.pdg, kStatusTarget, FourVector(nucleus.mass, 0, 0, 0));
  event.Add(hit.pdg, kStatusStruckNucleon, pBound, 1);
  event.Add(leptonPdg, kStatusFinal, lepton, 0);
  if (initial.remnantPdg != 0) event.Add(initial.remnantPdg, kStatusRemnant, initial.remnant, 1);
  return weight;
}

void DISKinematics::Finalize(const Nucleus& nucleus, Event& event, Random& rng) const {
  // Particles 0: beam, 2: struck nucleon, 3: lepton (see Generate).
  const FourVector k = event.particles[0].p;
  const FourVector pBound = event.particles[2].p;
  const FourVector lepton = event.particles[3].p;
  const FourVector X = k + pBound - lepton;
  const ThreeVector axis = (k - lepton).Vect();
  std::vector<DISHadronization::Hadron> hadrons;
  // Pion-multiplicity class: single pion with probability sigma_1pi / sigma_BY (theory/dis.md section 3),
  // otherwise two or more pions; within a class, the hadronization model's distribution.
  using PC = DISHadronization::PionClass;
  const int nuPdg = event.nuPdg, hitPdg = event.particles[2].pdg;
  const bool nu = nuPdg > 0, proton = hitPdg == pdg::kProton;
  const double W = X.M();
  const bool decomposed = opt_.isospinDecomposition || opt_.scaleI32 != 1.0 || opt_.scaleI12 != 1.0;
  PC cls = PC::kAny;
  SinglePionSplit split;
  if (decomposed) {
    const HitNucleon& hit = event.initial.hit[0];
    const double mIn = pdg::NucleonMass(hitPdg);
    const FourVector pOn(hit.OnShellEnergy(mIn), hit.p);
    const double lwSelf = ContractFor(nu, proton, k, lepton, pOn, mIn, W);
    const double lwPartner = NeedsPartner(nuPdg, hitPdg) ? ContractFor(nu, !proton, k, lepton, pOn, mIn, W) : 0;
    split = SinglePion(nuPdg, hitPdg, W, lwSelf, lwPartner);
    cls = rng.Uniform() < split.p1 ? PC::kSingle : PC::kMulti;
  }
  if (cls == PC::kSingle && opt_.isospinDecomposition) {
    std::vector<int> species;
    if (!NeedsPartner(nuPdg, hitPdg)) {
      species = nu ? std::vector<int>{pdg::kProton, pdg::kPiPlus} : std::vector<int>{pdg::kNeutron, pdg::kPiMinus};
    } else {
      const bool charged = rng.Uniform() < split.chargedShare;
      if (nu) species = charged ? std::vector<int>{pdg::kNeutron, pdg::kPiPlus} : std::vector<int>{pdg::kProton, pdg::kPiZero};
      else species = charged ? std::vector<int>{pdg::kProton, pdg::kPiMinus} : std::vector<int>{pdg::kNeutron, pdg::kPiZero};
    }
    had_.HadronizeSpecies(species, X, axis, rng, hadrons);
  } else {
    had_.Hadronize(nuPdg, hitPdg, X, axis, rng, hadrons, cls);
  }
  // Safety net: if the requested multiplicity class is closed (or not found), fall back to the
  // unrestricted hadronization.
  if (hadrons.empty() && cls != PC::kAny) had_.Hadronize(nuPdg, hitPdg, X, axis, rng, hadrons);
  if (hadrons.empty()) throw std::runtime_error("DISKinematics: hadronization failed");
  for (const auto& h : hadrons) event.Add(h.pdg, kStatusPreFSI, h.p, 2);
}

}  // namespace nuflow
