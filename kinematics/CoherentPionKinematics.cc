#include "kinematics/CoherentPionKinematics.h"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "core/Constants.h"
#include "core/PDG.h"

namespace nuflow {

using namespace constants;

namespace {

constexpr double kMbToGeV2 = 1.0 / 0.38937937;  // 1 mb = 2.568 GeV^-2
constexpr double kFmToInvGeV = 1.0 / 0.1973269804;

// BS09 Fig. 2 table: dsigma(pi 12C -> pi 12C)/dt = A1 exp(-b1 |t|) from the phase-shift analyses of the
// pi+- 12C elastic data (Dumbrajs et al., PRC 29, 581 (1984); Schlaile, PRC 55, 2584 (1997)).
struct ElasticNode {
  double T, A1, b1;  // GeV, mb/GeV^2, GeV^-2
};
constexpr ElasticNode kCarbon[] = {
    {0.076, 11600, 116.0}, {0.080, 14700, 109.0}, {0.100, 18300, 89.8}, {0.148, 21300, 91.0},
    {0.162, 22400, 89.2},  {0.226, 16400, 80.8},  {0.486, 5730, 54.6},  {0.584, 4610, 55.2},
    {0.662, 4570, 58.4},   {0.776, 4930, 60.5},   {0.870, 5140, 62.2},
};
constexpr int kNCarbon = sizeof(kCarbon) / sizeof(kCarbon[0]);

void CarbonTable(double T, double& A1, double& b1) {
  if (T <= kCarbon[0].T) {  // below the data: elastic strength vanishing linearly at threshold
    A1 = kCarbon[0].A1 * std::max(0.0, T) / kCarbon[0].T;
    b1 = kCarbon[0].b1;
    return;
  }
  if (T >= kCarbon[kNCarbon - 1].T) {
    A1 = kCarbon[kNCarbon - 1].A1;
    b1 = kCarbon[kNCarbon - 1].b1;
    return;
  }
  int i = 0;
  while (T > kCarbon[i + 1].T) ++i;
  const double f = (T - kCarbon[i].T) / (kCarbon[i + 1].T - kCarbon[i].T);
  A1 = kCarbon[i].A1 + f * (kCarbon[i + 1].A1 - kCarbon[i].A1);
  b1 = kCarbon[i].b1 + f * (kCarbon[i + 1].b1 - kCarbon[i].b1);
}

double PionLabS(double T) {
  const double E = T + kMassPionCharged;
  return kMassProton * kMassProton + kMassPionCharged * kMassPionCharged + 2 * kMassProton * E;
}

}  // namespace

double PionNucleusElastic::PiNTotal(double T) {
  // PDG (Review of Particle Physics, "Plots of cross sections and related quantities") Regge + ln^2 s fit
  // sigma(pi-+ p) = Z + B ln^2(s/sM) + Y1 (s1/s)^eta1 +- Y2 (s1/s)^eta2; the isoscalar average drops Y2.
  constexpr double Z = 20.86, Y1 = 19.24, eta1 = 0.4473;  // mb
  constexpr double M = 2.1206;                            // GeV
  const double B = kPi * 0.38937937 / (M * M);            // = 0.2720 mb
  const double sM = std::pow(kMassPionCharged + kMassProton + M, 2);
  const double s = PionLabS(T);
  const double l = std::log(s / sM);
  return Z + B * l * l + Y1 * std::pow(1.0 / s, eta1);
}

double PionNucleusElastic::PiNElastic(double T) {
  // Forward diffraction peak, sigma_el = sigma_tot^2 / (16 pi B_el) (real part neglected, as in BS09 Eq. 11),
  // with the effective slope B_el = B0 + 2 alpha' ln(s/s0), B0 = 5.5 GeV^-2, alpha' = 0.25 GeV^-2, s0 = 1 GeV^2
  // (pi p data at 5-200 GeV/c; Lebiedowicz & Szczurek, PRD 81, 036003 (2010), arXiv:0912.0190).
  const double s = PionLabS(T);
  const double Bel = 5.5 + 2 * 0.25 * std::log(std::max(s, 1.0));
  const double st = PiNTotal(T) * kMbToGeV2;
  return st * st / (16 * kPi * Bel) / kMbToGeV2;
}

void PionNucleusElastic::ReinSehgal(double T, int A, double& A1, double& b1) const {
  const double R0 = p_.R0 * kFmToInvGeV;                    // GeV^-1
  const double a13 = std::cbrt(static_cast<double>(A));
  b1 = R0 * R0 * a13 * a13 / 3.0;                           // BS09 Eq. 12
  const double sTot = PiNTotal(T) * kMbToGeV2;              // GeV^-2
  const double sInel = std::max(0.0, PiNTotal(T) - PiNElastic(T)) * kMbToGeV2;
  const double Fabs = std::exp(-9.0 * a13 / (16 * kPi * R0 * R0) * sInel);  // BS09 Eq. 13
  A1 = double(A) * A * sTot * sTot / (16 * kPi) * Fabs / kMbToGeV2;          // BS09 Eqs. 10-11 -> mb/GeV^2
}

void PionNucleusElastic::Coefficients(double T, int A, double& A1, double& b1) const {
  // 12C data fit; other nuclei: sigma_el ~ A^(2/3) (BS09 section 5), slope ~ R^2 ~ A^(2/3).
  auto table = [&](double TT, double& a, double& b) {
    CarbonTable(TT, a, b);
    const double r = std::cbrt(A / 12.0);
    a *= r * r * r * r;
    b *= r * r;
  };
  if (T <= p_.tableMaxT) return table(T, A1, b1);
  if (T >= p_.rsMinT) return ReinSehgal(T, A, A1, b1);
  double aT, bT, aR, bR;
  table(p_.tableMaxT, aT, bT);
  ReinSehgal(p_.rsMinT, A, aR, bR);
  const double f = (T - p_.tableMaxT) / (p_.rsMinT - p_.tableMaxT);
  A1 = aT + f * (aR - aT);
  b1 = bT + f * (bR - bT);
}

double PionNucleusElastic::DSigmaDt(double T, double absT, int A) const {
  double A1, b1;
  Coefficients(T, A, A1, b1);
  return A1 * std::exp(-b1 * absT);
}

double PionNucleusElastic::Sigma(double T, int A) const {
  double A1, b1;
  Coefficients(T, A, A1, b1);
  return A1 / b1;
}

std::string CoherentPionKinematics::Name() const {
  std::ostringstream s;
  s.precision(6);
  s << "CCCOH[BS09,fpi=" << opt_.fPi << ",mA=" << opt_.mA << ",lepmass=" << opt_.leptonMassCorrection
    << ",R0=" << elastic_.Params().R0 << ",T=" << elastic_.Params().tableMaxT << "-" << elastic_.Params().rsMinT
    << ",slope=" << opt_.proposalSlope << "]";
  return s.str();
}

double CoherentPionKinematics::DSigma(int nuPdg, double enu, double Q2, double nu, double absT, double Tpi,
                                      int A) const {
  const double mLep = pdg::LeptonMass(pdg::ChargedLeptonFor(nuPdg));
  const double mpi = kMassPionCharged;
  const double y = nu / enu;
  if (y <= 0 || y >= 1) return 0;
  const double ePrime = enu - nu;
  const double q3 = std::sqrt(nu * nu + Q2);
  const double u = (enu + ePrime + q3) / (2 * enu), v = (enu + ePrime - q3) / (2 * enu);
  if (v <= 0) return 0;
  const double GA = opt_.mA * opt_.mA / (opt_.mA * opt_.mA + Q2);
  double C = GA * GA;
  if (opt_.leptonMassCorrection && mLep > 0) {
    const double Q2min = mLep * mLep * y / (1 - y);
    const double d = Q2 + mpi * mpi;
    C = std::pow(GA - 0.5 * Q2min / d, 2) + 0.25 * y * std::max(0.0, Q2 - Q2min) * Q2min / (d * d);
  }
  const double dsdt = elastic_.DSigmaDt(Tpi, absT, A) * kMbToGeV2;  // GeV^-4
  return kGF * kGF * kVud * kVud * opt_.fPi * opt_.fPi / (2 * kPi * kPi) * (enu / q3) * u * v * C * dsdt;
}

double CoherentPionKinematics::Generate(int nuPdg, double enu, const Nucleus& nucleus, const InitialState& initial,
                                        const NuclearModel&, Event& event, Random& rng) const {
  const int leptonPdg = pdg::ChargedLeptonFor(nuPdg);
  const double mLep = pdg::LeptonMass(leptonPdg);
  const double mpi = kMassPionCharged;
  const int pionPdg = nuPdg > 0 ? pdg::kPiPlus : pdg::kPiMinus;
  const double MA = nucleus.mass;
  const int A = nucleus.A;

  // ---- nu = E - E': log-uniform in [m_pi, E - m_l].
  const double nuLo = mpi, nuHi = enu - mLep;
  if (nuHi <= nuLo * 1.0001) return 0;
  const double lr = std::log(nuHi / nuLo);
  const double nu = nuLo * std::exp(lr * rng.Uniform());
  const double gNu = 1.0 / (nu * lr);

  // ---- Q2 within the lepton kinematic limits; proposal ~ 1/(Q2 + mA^2)^2 (the G_A^2 shape).
  const double eP = enu - nu;
  const double pP = std::sqrt(std::max(0.0, eP * eP - mLep * mLep));
  if (pP <= 0) return 0;
  const double Q2lo = std::max(0.0, 2 * enu * (eP - pP) - mLep * mLep);
  const double Q2hi = 2 * enu * (eP + pP) - mLep * mLep;
  const double L2 = opt_.mA * opt_.mA;
  const double a = 1 / (Q2lo + L2), b = 1 / (Q2hi + L2);
  const double Q2 = 1 / (a - rng.Uniform() * (a - b)) - L2;
  const double gQ2 = 1 / ((Q2 + L2) * (Q2 + L2) * (a - b));

  // ---- |t|: exponential proposal above (0.9 x) the recoil-free collinear minimum (|q| - p_pi)^2.
  const double q3 = std::sqrt(nu * nu + Q2);
  const double p0 = std::sqrt(std::max(0.0, nu * nu - mpi * mpi));
  const double tMin0 = std::pow((Q2 + mpi * mpi) / (q3 + p0), 2);
  const double tLo = 0.9 * tMin0;
  const double R0 = elastic_.Params().R0 / 0.1973269804;
  const double a13 = std::cbrt(static_cast<double>(A));
  const double slope = opt_.proposalSlope * R0 * R0 * a13 * a13 / 3.0;
  const double absT = tLo - std::log(rng.UniformOpen()) / slope;
  const double gT = slope * std::exp(-slope * (absT - tLo));

  // ---- Pion: E_pi = nu - |t|/(2 M_A) (elastic recoil of the ground-state nucleus), angle to q from t.
  const double ePi = nu - absT / (2 * MA);
  if (ePi <= mpi) return 0;
  const double pPi = std::sqrt(ePi * ePi - mpi * mpi);
  const double cosPq = (-absT + Q2 - mpi * mpi + 2 * nu * ePi) / (2 * q3 * pPi);
  if (std::abs(cosPq) > 1) return 0;

  const double ds = DSigma(nuPdg, enu, Q2, nu, absT, ePi - mpi, A);  // per dQ2 dy d|t|
  if (ds <= 0) return 0;
  const double weight = ds / enu / (gNu * gQ2 * gT) * kGeV2To1e38Cm2;  // dy = dnu / E

  // ---- Four-vectors (neutrino along z).
  const double cosL = std::clamp((2 * enu * eP - mLep * mLep - Q2) / (2 * enu * pP), -1.0, 1.0);
  const double sinL = std::sqrt(std::max(0.0, 1 - cosL * cosL));
  const double phiL = 2 * kPi * rng.Uniform();
  const FourVector k(enu, 0, 0, enu);
  const FourVector lepton(eP, pP * sinL * std::cos(phiL), pP * sinL * std::sin(phiL), pP * cosL);
  const FourVector q = k - lepton;
  const ThreeVector zq = q.Vect().Unit();
  const ThreeVector ref = std::abs(zq.x) < 0.9 ? ThreeVector{1, 0, 0} : ThreeVector{0, 1, 0};
  const ThreeVector xq = ref.Cross(zq).Unit();
  const ThreeVector yq = zq.Cross(xq);
  const double sinPq = std::sqrt(std::max(0.0, 1 - cosPq * cosPq));
  const double phiP = 2 * kPi * rng.Uniform();
  const ThreeVector dir = zq * cosPq + xq * (sinPq * std::cos(phiP)) + yq * (sinPq * std::sin(phiP));
  const FourVector pion(ePi, dir * pPi);
  const FourVector target(MA, 0, 0, 0);
  const FourVector recoil = q + target - pion;

  event.nuPdg = nuPdg;
  event.targetPdg = nucleus.pdg;
  event.enu = enu;
  event.Q2 = Q2;
  event.q0 = nu;
  event.q3 = q3;
  const double M = kMassNucleon;
  event.W = std::sqrt(std::max(0.0, M * M + 2 * M * nu - Q2));  // free-nucleon-at-rest convention
  event.initial = initial;
  event.particles.clear();
  event.Add(nuPdg, kStatusBeam, k);
  event.Add(nucleus.pdg, kStatusTarget, target);
  event.Add(leptonPdg, kStatusFinal, lepton, 0);
  // The coherent pion leaves the nucleus directly: pion distortion and absorption are part of the
  // measured elastic pi-A amplitude, so it is not handed to the cascade (status 1, not 22).
  event.Add(pionPdg, kStatusFinal, pion, 1);
  // Recoiling ground-state nucleus (|t| = -(P_A' - P_A)^2); written as the remnant so that it is never a
  // status-1 hadron (the cascade re-emits it as the status-25 remnant with the same four-momentum).
  event.Add(nucleus.pdg, kStatusRemnant, recoil, 1);
  return weight;
}

}  // namespace nuflow
