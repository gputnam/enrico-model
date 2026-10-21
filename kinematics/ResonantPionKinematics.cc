#include "kinematics/ResonantPionKinematics.h"

#include <algorithm>
#include <cmath>
#include <sstream>

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

ThreeVector IsotropicDirection(Random& rng) {
  const double c = 2 * rng.Uniform() - 1;
  const double s = std::sqrt(std::max(0.0, 1 - c * c));
  const double phi = 2 * kPi * rng.Uniform();
  return {s * std::cos(phi), s * std::sin(phi), c};
}

struct FinalState {
  int nucleon, pion;
  PionIsospin iso;
};

constexpr double kBWWidth = 0.12;  // width of the W proposal Breit-Wigner [GeV] (proposal only)

}  // namespace

std::string ResonantPionKinematics::Name() const {
  std::ostringstream s;
  s.precision(6);
  s << "CCRES[" << current_.Name() << ",Wmax=" << opt_.wMax << "]";
  return s.str();
}

double ResonantPionKinematics::Generate(int nuPdg, double enu, const Nucleus& nucleus, const InitialState& initial,
                                        const NuclearModel& nuclearModel, Event& event, Random& rng) const {
  const HitNucleon& hit = initial.hit[0];
  const bool nu = nuPdg > 0;
  const int helicity = nu ? +1 : -1;
  const int leptonPdg = pdg::ChargedLeptonFor(nuPdg);
  const double mLep = pdg::LeptonMass(leptonPdg);
  const double mIn = pdg::NucleonMass(hit.pdg);

  // Final states allowed by charge conservation, with the isospin amplitude that describes them.
  FinalState fs[2];
  int nFS = 0;
  if (nu) {
    if (hit.pdg == pdg::kProton) {
      fs[nFS++] = {pdg::kProton, pdg::kPiPlus, PionIsospin::PPiPlus};
    } else {
      fs[nFS++] = {pdg::kNeutron, pdg::kPiPlus, PionIsospin::NPiPlus};
      fs[nFS++] = {pdg::kProton, pdg::kPiZero, PionIsospin::PPiZero};
    }
  } else {
    if (hit.pdg == pdg::kNeutron) {
      fs[nFS++] = {pdg::kNeutron, pdg::kPiMinus, PionIsospin::PPiPlus};
    } else {
      fs[nFS++] = {pdg::kProton, pdg::kPiMinus, PionIsospin::NPiPlus};
      fs[nFS++] = {pdg::kNeutron, pdg::kPiZero, PionIsospin::PPiZero};
    }
  }
  const FinalState& f = fs[nFS == 1 ? 0 : (rng.Uniform() < 0.5 ? 0 : 1)];
  const double mN = pdg::NucleonMass(f.nucleon);
  const double mPi = pdg::PionMass(f.pion);

  const FourVector k(enu, 0, 0, enu);
  const FourVector pBound = hit.Momentum();
  const FourVector Ptot = k + pBound;
  const double s = Ptot.M2();
  if (Ptot.e <= 0 || s <= 0) return 0;
  const double sqrtS = std::sqrt(s);
  const double wLo = mN + mPi;
  const double wHi = std::min(opt_.wMax, sqrtS - mLep);
  if (wHi <= wLo) return 0;

  // --- W proposal: Breit-Wigner (Delta) + flat mixture on [wLo, wHi].
  const double mD = current_.Parameters().MDelta;
  const double aLo = std::atan(2 * (wLo - mD) / kBWWidth);
  const double aHi = std::atan(2 * (wHi - mD) / kBWWidth);
  double W;
  if (rng.Uniform() < opt_.bwFraction) {
    W = mD + 0.5 * kBWWidth * std::tan(aLo + rng.Uniform() * (aHi - aLo));
    W = std::clamp(W, wLo, wHi);
  } else {
    W = wLo + rng.Uniform() * (wHi - wLo);
  }
  const double xBW = 2 * (W - mD) / kBWWidth;
  const double gW = opt_.bwFraction * (2 / kBWWidth) / (1 + xBW * xBW) / (aHi - aLo) +
                    (1 - opt_.bwFraction) / (wHi - wLo);

  // --- Lepton in the CM frame of k + p_bound, against the hadronic system of mass W.
  const ThreeVector beta = Ptot.BoostVector();
  const FourVector kCM = k.Boosted(-beta);
  const double kStar = kCM.e;
  const double lamL = Kallen(s, mLep * mLep, W * W);
  if (lamL <= 0) return 0;
  const double pLepStar = std::sqrt(lamL) / (2 * sqrtS);
  const double eLepStar = std::sqrt(pLepStar * pLepStar + mLep * mLep);

  const double bCoef = 2 * kStar * pLepStar / opt_.proposalQ2Scale;
  double x, gx;
  if (bCoef < 1e-6) {
    x = 2 * rng.Uniform();
    gx = 0.5;
  } else {
    const double kappa = 2 * bCoef / (1 + 2 * bCoef);
    x = (1.0 / (1.0 - rng.Uniform() * kappa) - 1.0) / bCoef;
    const double d = 1 + bCoef * x;
    gx = bCoef / (kappa * d * d);
  }
  const double cosT = std::clamp(1.0 - x, -1.0, 1.0);
  const double sinT = std::sqrt(std::max(0.0, 1 - cosT * cosT));
  const double phi = 2 * kPi * rng.Uniform();
  const ThreeVector n = kCM.Vect().Unit();
  ThreeVector e1, e2;
  Perpendicular(n, e1, e2);
  const ThreeVector dir = n * cosT + e1 * (sinT * std::cos(phi)) + e2 * (sinT * std::sin(phi));
  const FourVector lepton = FourVector(eLepStar, dir * pLepStar).Boosted(beta);
  const FourVector X = Ptot - lepton;  // hadronic system, mass W

  // --- X -> N pi, isotropic proposal in the X rest frame.
  const double lamH = Kallen(W * W, mN * mN, mPi * mPi);
  if (lamH <= 0) return 0;
  const double pPiStar = std::sqrt(lamH) / (2 * W);
  const ThreeVector dPi = IsotropicDirection(rng);
  const ThreeVector betaX = X.BoostVector();
  const FourVector pion = FourVector(std::sqrt(pPiStar * pPiStar + mPi * mPi), dPi * pPiStar).Boosted(betaX);
  const FourVector nucleonOut = X - pion;

  if (nuclearModel.PauliBlocked(nucleus, f.nucleon, nucleonOut.Vect())) return 0;

  // --- Hadronic vertex with the on-shell initial nucleon (de Forest): q~ = p' + k_pi - p_on.
  const double eOn = hit.OnShellEnergy(mIn);
  const FourVector pOn(eOn, hit.p);
  const FourVector qTilde = nucleonOut + pion - pOn;
  if (-qTilde.M2() <= 0) return 0;

  PionCurrentPieces pieces;
  current_.Pieces(pOn, nucleonOut, pion, pieces);
  dirac::Mat S[4];
  current_.Combine(pieces, f.iso, S);
  const double LH = PionProductionCurrent::Contract(k, lepton, mLep, helicity, pOn, mIn, nucleonOut, mN, S);
  if (!(LH > 0)) return 0;

  // dsigma = <|M|^2> / (4 E_nu E_N) dPhi_3,  <|M|^2> = G_F^2 Vud^2 LH / 4,
  // dPhi_3 = [p_l*/(16 pi^2 sqrt s) dOmega_l*] [dW^2 / 2pi] [p_pi*/(16 pi^2 W) dOmega_pi*].
  const double msq = kGF * kGF * kVud * kVud * 0.25 * LH;
  const double phaseSpace = pLepStar / (16 * kPi * kPi * sqrtS) * (2 * W / (2 * kPi)) * pPiStar / (16 * kPi * kPi * W);
  const double invProposal = (2 * kPi / gx) * (4 * kPi) / gW * nFS;
  const double weight = msq / (4 * enu * eOn) * phaseSpace * invProposal * kGeV2To1e38Cm2;

  // --- Primary vertex.
  const FourVector q = k - lepton;
  event.nuPdg = nuPdg;
  event.targetPdg = nucleus.pdg;
  event.enu = enu;
  event.Q2 = -q.M2();
  event.q0 = q.e;
  event.q3 = q.P();
  event.W = W;
  event.initial = initial;
  auto& parts = event.particles;
  parts.clear();
  event.Add(nuPdg, kStatusBeam, k);
  event.Add(nucleus.pdg, kStatusTarget, FourVector(nucleus.mass, 0, 0, 0));
  const int iHit = event.Add(hit.pdg, kStatusStruckNucleon, pBound, 1);
  event.Add(leptonPdg, kStatusFinal, lepton, 0);
  event.Add(f.nucleon, kStatusPreFSI, nucleonOut, iHit);
  event.Add(f.pion, kStatusPreFSI, pion, iHit);
  if (initial.remnantPdg != 0) event.Add(initial.remnantPdg, kStatusRemnant, initial.remnant, 1);
  return weight;
}

}  // namespace nuflow
