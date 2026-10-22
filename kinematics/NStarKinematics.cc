#include "kinematics/NStarKinematics.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "core/Constants.h"
#include "core/Dirac.h"
#include "core/PDG.h"
#include "kinematics/DISHadronization.h"

namespace nuflow {

using namespace constants;
using dirac::cplx;
using dirac::Comp;
using dirac::Gammas;
using dirac::Mat;
using dirac::Metric;

namespace {

constexpr int kEta = 221;
constexpr double kMassEta = 0.547862;  // PDG 2024 [GeV]

double Kallen(double a, double b, double c) { return a * a + b * b + c * c - 2 * (a * b + a * c + b * c); }

double PStar(double W, double m1, double m2) {
  if (W <= m1 + m2) return 0;
  return std::sqrt(std::max(0.0, Kallen(W * W, m1 * m1, m2 * m2))) / (2 * W);
}

void Perpendicular(const ThreeVector& n, ThreeVector& e1, ThreeVector& e2) {
  const ThreeVector a = std::abs(n.x) < 0.9 ? ThreeVector{1, 0, 0} : ThreeVector{0, 1, 0};
  e1 = n.Cross(a).Unit();
  e2 = n.Cross(e1);
}

// Two-body phase-space factor with centrifugal barrier: p^(2L+1) / W.
double Rho2(double W, double m1, double m2, int L) {
  const double p = PStar(W, m1, m2);
  return p > 0 ? std::pow(p, 2 * L + 1) / W : 0.0;
}

// Three-body (N pi pi) phase-space volume up to a constant: (1/W) int dm p(W -> N + m) p(m -> pi pi).
double Rho3(double W) {
  const double mN = kMassNucleon, mpi = kMassPionAverage;
  const double lo = 2 * mpi, hi = W - mN;
  if (hi <= lo) return 0;
  constexpr int n = 48;  // Simpson
  const double h = (hi - lo) / n;
  double s = 0;
  for (int i = 0; i <= n; ++i) {
    const double m = lo + i * h;
    const double f = PStar(W, mN, m) * PStar(m, mpi, mpi);
    s += f * (i == 0 || i == n ? 1 : (i % 2 ? 4 : 2));
  }
  return s * h / 3 / W;
}

constexpr double kProposalFlat = 0.25;  // fraction of flat W proposals

}  // namespace

NStarResonance NStarKinematics::Lookup(const std::string& name) {
  // LPP Table I (mass, total width, elasticity x_pi). The inelastic remainder: N(1440), N(1520) -> N pi pi;
  // N(1535) -> N eta 0.45 and N pi pi 0.15 (PDG 2024 ranges: N eta 30-55%, N pi pi 4-31%).
  if (name == "P11_1440") return {name, 1, +1, 1.440, 0.350, 1, 0.60, 0.00, 0.40};
  if (name == "D13_1520") return {name, 3, -1, 1.520, 0.125, 2, 0.50, 0.00, 0.50};
  if (name == "S11_1535") return {name, 1, -1, 1.535, 0.150, 0, 0.40, 0.45, 0.15};
  throw std::runtime_error("CCNSTAR: unknown resonance '" + name + "' (P11_1440, D13_1520, S11_1535)");
}

NStarKinematics::NStarKinematics(Options opt) : opt_(std::move(opt)) {
  for (const auto& n : opt_.resonances) res_.push_back(Lookup(n));
  if (res_.empty() || res_.size() > 8) throw std::runtime_error("CCNSTAR: 1 to 8 resonances required");
}

std::string NStarKinematics::Name() const {
  std::ostringstream s;
  s.precision(6);
  s << "CCNSTAR[LPP2006,";
  for (const auto& r : res_) s << r.name << ",";
  s << "MV=" << opt_.MV << ",MA=" << opt_.MA << ",Wmax=" << opt_.wMax;
  if (opt_.vectorScale != 1.0 || opt_.axialScale != 1.0) s << ",scaleVA=" << opt_.vectorScale << "/" << opt_.axialScale;
  s << "]";
  return s.str();
}

NStarFormFactors NStarKinematics::FormFactors(const NStarResonance& r, double Q2) const {
  const double m = kMassNucleon, mpi2 = kMassPionAverage * kMassPionAverage;
  const double MV2 = opt_.MV * opt_.MV, MA2 = opt_.MA * opt_.MA;
  const double DV = (1 + Q2 / MV2) * (1 + Q2 / MV2), DA = (1 + Q2 / MA2) * (1 + Q2 / MA2);
  const double lnQ = std::log(1 + Q2);  // ln(1 + Q2 / 1 GeV^2)
  const double axialShape = 1 / DA / (1 + Q2 / (3 * MA2));
  NStarFormFactors f;
  if (r.name == "P11_1440") {
    // LPP Eq. (4.13): proton e.m. form factors; isoscalar neglected, g^V = -2 g^(p).
    const double g1p = 2.3 / DV / (1 + Q2 / (4.3 * MV2));
    const double g2p = -0.76 / DV * (1 - 2.8 * lnQ);
    f.g1V = -2 * g1p;
    f.g2V = -2 * g2p;
    f.g1A = -0.51 * axialShape;
    f.g3A = (r.M + m) * m / (Q2 + mpi2) * f.g1A;
  } else if (r.name == "S11_1535") {
    // LPP Eq. (4.18).
    const double g1p = 2.0 / DV / (1 + Q2 / (1.2 * MV2)) * (1 + 7.2 * lnQ);
    const double g2p = 0.84 / DV * (1 + 0.11 * lnQ);
    f.g1V = -2 * g1p;
    f.g2V = -2 * g2p;
    f.g1A = -0.21 * axialShape;
    f.g3A = (r.M - m) * m / (Q2 + mpi2) * f.g1A;
  } else if (r.name == "D13_1520") {
    // LPP Eq. (4.5): C_i^V = C_i^n - C_i^p; C5A(0) = -2.1 (PCAC), C6A pion pole, C3A = C4A = 0.
    const double s = 1 / DV / (1 + Q2 / (8.9 * MV2));
    f.C3V = (-1.13 - 2.95) * s;
    f.C4V = (0.46 + 1.05) * s;
    f.C5V = (-0.17 + 0.48) / DV;
    f.C5A = -2.1 * axialShape;
    f.C6A = m * m * f.C5A / (mpi2 + Q2);
  }
  f.g1V *= opt_.vectorScale;
  f.g2V *= opt_.vectorScale;
  f.C3V *= opt_.vectorScale;
  f.C4V *= opt_.vectorScale;
  f.C5V *= opt_.vectorScale;
  f.g1A *= opt_.axialScale;
  f.g3A *= opt_.axialScale;
  f.C5A *= opt_.axialScale;
  f.C6A *= opt_.axialScale;
  return f;
}

void NStarKinematics::PartialWidths(const NStarResonance& r, double W, double& gPi, double& gEta,
                                    double& gPiPi) const {
  const double mN = kMassNucleon, mpi = kMassPionAverage;
  const double rPi0 = Rho2(r.M, mN, mpi, r.lPi);
  gPi = r.xPi > 0 ? r.Gamma0 * r.xPi * Rho2(W, mN, mpi, r.lPi) / rPi0 : 0;
  gEta = r.xEta > 0 ? r.Gamma0 * r.xEta * Rho2(W, mN, kMassEta, 0) / Rho2(r.M, mN, kMassEta, 0) : 0;
  gPiPi = r.xPiPi > 0 ? r.Gamma0 * r.xPiPi * Rho3(W) / Rho3(r.M) : 0;
}

double NStarKinematics::Spectral(const NStarResonance& r, double W) const {
  double a, b, c;
  PartialWidths(r, W, a, b, c);
  const double G = a + b + c;
  const double d = W * W - r.M * r.M;
  return r.M * G / kPi / (d * d + r.M * r.M * G * G);
}

double NStarKinematics::ContractLH(const NStarResonance& r, const FourVector& k, const FourVector& kp, double mLep,
                                   int helicity, const FourVector& p, double mIn, const FourVector& pR) const {
  const Gammas& G = Gammas::Get();
  const double m = kMassNucleon, MR = r.M;
  const FourVector q = pR - p;
  const double Q2 = -q.M2();
  const NStarFormFactors f = FormFactors(r, Q2);
  const Mat qs = G.Slash(q);
  const Mat Pi = G.Slash(p) + G.one * mIn;

  // Lepton tensor (as PionProductionCurrent::Contract).
  static const Mat PL = Gammas::Get().one - Gammas::Get().g5;
  const Mat Kin = helicity > 0 ? G.Slash(k) : G.Slash(kp) - G.one * mLep;
  const Mat Kout = helicity > 0 ? G.Slash(kp) + G.one * mLep : G.Slash(k);
  cplx L[4][4];
  {
    Mat U[4], gPL[4];
    for (int mu = 0; mu < 4; ++mu) {
      gPL[mu] = G.g[mu] * PL;
      U[mu] = Kout * gPL[mu] * Kin;
    }
    for (int mu = 0; mu < 4; ++mu)
      for (int nu = 0; nu < 4; ++nu) L[mu][nu] = dirac::TraceProduct(U[mu], gPL[nu]);
  }

  cplx H[4][4];
  if (r.twoJ == 1) {
    // <R| J^nu |N> = ubar(pR) Gamma^nu u(p) (LPP Eqs. 4.10, 4.15), mu_R = m_N + M_R.
    const double muR = m + MR;
    Mat S[4];
    for (int nu = 0; nu < 4; ++nu) {
      Mat V = G.g[nu] * (f.g1V * Q2 / (muR * muR));
      V.AddScaled(qs, f.g1V * Comp(q, nu) / (muR * muR));
      // i sigma^{nu rho} q_rho = -(1/2)(gamma^nu qslash - qslash gamma^nu)
      V.AddScaled(G.g[nu] * qs - qs * G.g[nu], -0.5 * f.g2V / muR);
      if (r.parity > 0) {
        S[nu] = V;
        S[nu].AddScaled(G.g5g[nu], -f.g1A);
        S[nu].AddScaled(G.g5, -f.g3A * Comp(q, nu) / m);
      } else {
        S[nu] = V * G.g5;
        S[nu].AddScaled(G.g[nu], -f.g1A);
        S[nu].AddScaled(G.one, -f.g3A * Comp(q, nu) / m);
      }
    }
    const Mat PR = G.Slash(pR) + G.one * MR;
    Mat A[4], B[4];
    for (int mu = 0; mu < 4; ++mu) {
      A[mu] = PR * S[mu] * Pi;
      B[mu] = G.Bar(S[mu]);
    }
    for (int mu = 0; mu < 4; ++mu)
      for (int nu = 0; nu < 4; ++nu) H[mu][nu] = dirac::TraceProduct(A[mu], B[nu]);
  } else {
    // <R| J^nu |N> = psibar_lambda(pR) d^{lambda nu} u(p) (LPP Eq. 4.2), C3A = C4A = 0:
    //   V^{lambda nu} = g^{lambda nu}[C3V qslash/m + C4V pR.q/m^2 + C5V p.q/m^2] - q^lambda[C3V gamma^nu/m + C4V pR^nu/m^2 + C5V p^nu/m^2]
    //   A^{lambda nu} = g^{lambda nu} C5A + q^lambda q^nu C6A/m^2
    //   d = V + A gamma5 (negative parity, D13);  d = V gamma5 + A (positive parity, P33-like).
    const double m2 = m * m;
    const double pRq = pR.Dot(q), pq = p.Dot(q);
    Mat D[4][4];  // D[lambda][nu]
    for (int la = 0; la < 4; ++la)
      for (int nu = 0; nu < 4; ++nu) {
        Mat V;
        const double g = (la == nu) ? Metric(la) : 0.0;
        if (g != 0) {
          V.AddScaled(qs, g * f.C3V / m);
          V.AddScaled(G.one, g * (f.C4V * pRq + f.C5V * pq) / m2);
        }
        V.AddScaled(G.g[nu], -Comp(q, la) * f.C3V / m);
        V.AddScaled(G.one, -Comp(q, la) * (f.C4V * Comp(pR, nu) + f.C5V * Comp(p, nu)) / m2);
        const double a = g * f.C5A + Comp(q, la) * Comp(q, nu) * f.C6A / m2;
        if (r.parity < 0) {
          D[la][nu] = V;
          D[la][nu].AddScaled(G.g5, a);
        } else {
          D[la][nu] = V * G.g5;
          D[la][nu].AddScaled(G.one, a);
        }
      }
    // Rarita-Schwinger spin sum P^{rho lambda}(pR) = -(pRslash + M)[g - gamma gamma/3 - 2 pR pR/(3M^2) + (pR^rho gamma^lambda - pR^lambda gamma^rho)/(3M)]
    const Mat pre = (G.Slash(pR) + G.one * MR) * (-1.0);
    Mat Y[4][4];  // Y[rho][mu] = sum_lambda g_ll P^{rho lambda} D^{lambda mu} (pslash + m)
    Mat P[4][4];
    for (int rho = 0; rho < 4; ++rho)
      for (int la = 0; la < 4; ++la) {
        Mat B = G.one * ((rho == la ? Metric(rho) : 0.0) - 2 * Comp(pR, rho) * Comp(pR, la) / (3 * MR * MR));
        B.AddScaled(G.g[rho] * G.g[la], -1.0 / 3);
        B.AddScaled(G.g[la], Comp(pR, rho) / (3 * MR));
        B.AddScaled(G.g[rho], -Comp(pR, la) / (3 * MR));
        P[rho][la] = pre * B;
      }
    for (int rho = 0; rho < 4; ++rho)
      for (int mu = 0; mu < 4; ++mu) {
        Mat s;
        for (int la = 0; la < 4; ++la) s.AddScaled(P[rho][la] * D[la][mu], Metric(la));
        Y[rho][mu] = s * Pi;
      }
    Mat Db[4][4];
    for (int rho = 0; rho < 4; ++rho)
      for (int nu = 0; nu < 4; ++nu) Db[rho][nu] = G.Bar(D[rho][nu]);
    for (int mu = 0; mu < 4; ++mu)
      for (int nu = 0; nu < 4; ++nu) {
        cplx h = 0;
        for (int rho = 0; rho < 4; ++rho) h += Metric(rho) * dirac::TraceProduct(Y[rho][mu], Db[rho][nu]);
        H[mu][nu] = h;
      }
  }

  double total = 0;
  for (int mu = 0; mu < 4; ++mu)
    for (int nu = 0; nu < 4; ++nu) total += (L[mu][nu] * H[mu][nu]).real() * Metric(mu) * Metric(nu);
  return total;
}

double NStarKinematics::Generate(int nuPdg, double enu, const Nucleus& nucleus, const InitialState& initial,
                                 const NuclearModel& nuclearModel, Event& event, Random& rng) const {
  const HitNucleon& hit = initial.hit[0];
  const bool nu = nuPdg > 0;
  // Charged-current excitation of an I = 1/2 resonance: W+ n -> R+, W- p -> R0 only.
  if (nu != (hit.pdg == pdg::kNeutron)) return 0;
  const int helicity = nu ? +1 : -1;
  const int leptonPdg = pdg::ChargedLeptonFor(nuPdg);
  const double mLep = pdg::LeptonMass(leptonPdg);
  const double mIn = pdg::NucleonMass(hit.pdg);

  const FourVector k(enu, 0, 0, enu);
  const FourVector pBound = hit.Momentum();
  const FourVector Ptot = k + pBound;
  const double s = Ptot.M2();
  if (Ptot.e <= 0 || s <= 0) return 0;
  const double sqrtS = std::sqrt(s);
  const double wLo = (nu ? kMassProton : kMassNeutron) + kMassPionNeutral;
  const double wHi = std::min(opt_.wMax, sqrtS - mLep);
  if (wHi <= wLo) return 0;

  // --- W proposal: flat + one truncated Cauchy per resonance.
  const int nR = static_cast<int>(res_.size());
  double W;
  if (rng.Uniform() < kProposalFlat) {
    W = wLo + rng.Uniform() * (wHi - wLo);
  } else {
    const auto& r = res_[std::min(nR - 1, static_cast<int>(rng.Uniform() * nR))];
    const double h = 0.5 * r.Gamma0;
    const double aLo = std::atan((wLo - r.M) / h), aHi = std::atan((wHi - r.M) / h);
    W = std::clamp(r.M + h * std::tan(aLo + rng.Uniform() * (aHi - aLo)), wLo, wHi);
  }
  double gW = kProposalFlat / (wHi - wLo);
  for (const auto& r : res_) {
    const double h = 0.5 * r.Gamma0;
    const double aLo = std::atan((wLo - r.M) / h), aHi = std::atan((wHi - r.M) / h);
    const double x = (W - r.M) / h;
    gW += (1 - kProposalFlat) / nR / (h * (1 + x * x) * (aHi - aLo));
  }

  // --- Lepton in the CM frame of k + p_bound (as CCRES).
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
  const FourVector X = Ptot - lepton;  // resonance, mass W

  // --- Hadronic vertex with the on-shell initial nucleon (de Forest), q~ = X - p_on.
  const double eOn = hit.OnShellEnergy(mIn);
  const FourVector pOn(eOn, hit.p);
  if (-(X - pOn).M2() <= 0) return 0;
  double contrib[8];
  double sum = 0;
  for (int i = 0; i < nR; ++i) {
    const double lh = ContractLH(res_[i], k, lepton, mLep, helicity, pOn, mIn, X);
    contrib[i] = lh > 0 ? lh * Spectral(res_[i], W) : 0;
    sum += contrib[i];
  }
  if (!(sum > 0)) return 0;

  // dsigma = <|M|^2> / (4 E_nu E_N) dPhi_2 A(W^2) dW^2,  <|M|^2> = G_F^2 Vud^2 LH / 4,
  // dPhi_2 = p_l* / (16 pi^2 sqrt s) dOmega_l*.
  const double weight = kGF * kGF * kVud * kVud * 0.25 * sum / (4 * enu * eOn) * pLepStar / (16 * kPi * kPi * sqrtS) *
                        (2 * W) * (2 * kPi / gx) / gW * kGeV2To1e38Cm2;

  // --- Resonance and decay channel.
  double u = rng.Uniform() * sum;
  int iR = 0;
  while (iR < nR - 1 && u >= contrib[iR]) u -= contrib[iR++];
  const NStarResonance& r = res_[iR];
  double gPi, gEta, gPiPi;
  PartialWidths(r, W, gPi, gEta, gPiPi);
  const double g = rng.Uniform() * (gPi + gEta + gPiPi);
  // Charge states for R+ (nu n); R0 (nubar p) is the isospin mirror (p <-> n, pi+ <-> pi-).
  //  N pi:    n pi+ 2/3, p pi0 1/3.
  //  N eta:   p eta.
  //  N pi pi via pi Delta (I = 1/2): p pi+ pi- 5/9, p pi0 pi0 2/9, n pi+ pi0 2/9.
  std::vector<int> species;
  const double v = rng.Uniform();
  if (g < gPi) {
    species = v < 2.0 / 3 ? std::vector<int>{pdg::kNeutron, pdg::kPiPlus} : std::vector<int>{pdg::kProton, pdg::kPiZero};
  } else if (g < gPi + gEta) {
    species = {pdg::kProton, kEta};
  } else {
    if (v < 5.0 / 9) species = {pdg::kProton, pdg::kPiPlus, pdg::kPiMinus};
    else if (v < 7.0 / 9) species = {pdg::kProton, pdg::kPiZero, pdg::kPiZero};
    else species = {pdg::kNeutron, pdg::kPiPlus, pdg::kPiZero};
  }
  if (!nu)
    for (int& c : species) {
      if (c == pdg::kProton) c = pdg::kNeutron;
      else if (c == pdg::kNeutron) c = pdg::kProton;
      else if (c == pdg::kPiPlus) c = pdg::kPiMinus;
      else if (c == pdg::kPiMinus) c = pdg::kPiPlus;
    }
  std::vector<double> masses;
  for (int c : species)
    masses.push_back(c == pdg::kProton || c == pdg::kNeutron ? pdg::NucleonMass(c)
                     : c == kEta                            ? kMassEta
                                                            : pdg::PionMass(c));
  std::vector<FourVector> pd;
  if (!DISHadronization::PhaseSpace(W, masses, rng, pd)) return 0;  // within MeV of a charge-state threshold
  const ThreeVector betaX = X.BoostVector();
  for (auto& v4 : pd) v4 = v4.Boosted(betaX);
  if (nuclearModel.PauliBlocked(nucleus, species[0], pd[0].Vect())) return 0;

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
  for (size_t i = 0; i < species.size(); ++i) event.Add(species[i], kStatusPreFSI, pd[i], iHit);
  if (initial.remnantPdg != 0) event.Add(initial.remnantPdg, kStatusRemnant, initial.remnant, 1);
  return weight;
}

}  // namespace nuflow
