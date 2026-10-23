#include "kinematics/ValenciaRPA.h"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "core/Constants.h"

namespace nuflow {

using constants::kPi;

namespace {

constexpr double kHbarC = 0.1973269804;  // GeV fm
constexpr double kMPi = 0.13957;         // charged pion mass used in V_l (Nieves et al. 2004) [GeV]
constexpr double kMN = 0.938919;         // average nucleon mass [GeV]

// Integral  int_0^kF p ln|(a + p)/(a - p)| dp.
double LindhardF(double a, double kF) {
  // The logarithmic term vanishes (x ln x -> 0) at a = +-kF.
  if (std::abs(std::abs(a) - kF) < 1e-12 * kF) return a * kF;
  return a * kF + 0.5 * (kF * kF - a * a) * std::log(std::abs((a + kF) / (a - kF)));
}

enum Block { k00 = 0, k0L, kLL, kTT, k0T, kLT, kNBlocks };

inline Block BlockOf(int m, int n) {
  const bool t0 = m == 0, t1 = n == 0;
  const bool l0 = m == 3, l1 = n == 3;
  if (t0 && t1) return k00;
  if ((t0 && l1) || (l0 && t1)) return k0L;
  if (l0 && l1) return kLL;
  if (t0 || t1) return k0T;
  if (l0 || l1) return kLT;
  return kTT;
}

// Components (e, x, y, z) of a four-vector in the frame whose z axis is ez.
inline void Components(const FourVector& v, const ThreeVector& ex, const ThreeVector& ey, const ThreeVector& ez,
                       double out[4]) {
  const ThreeVector p = v.Vect();
  out[0] = v.e;
  out[1] = p.Dot(ex);
  out[2] = p.Dot(ey);
  out[3] = p.Dot(ez);
}

}  // namespace

std::string ValenciaRPA::Describe() const {
  std::ostringstream s;
  s << "RPA(g'=" << par_.gPrime << ",f2=" << par_.f2Over4Pi << ",Lpi=" << par_.lambdaPi << ",Crho=" << par_.cRho
    << ",Lrho=" << par_.lambdaRho << ",f*/f=" << par_.fStarOverF << ",GD=" << par_.gammaDelta << ",Gph=" << par_.gammaPH << ",c0=" << par_.c0 << ",f'=" << par_.f0pIn << "/"
    << par_.f0pEx << (par_.delta ? ",Dh" : "") << ",q=" << par_.qFull << "-" << par_.qZero << ")";
  return s.str();
}

std::complex<double> ValenciaRPA::LindhardN(double q0, double q, double kF, double M) {
  if (kF <= 0 || q <= 0) return {0.0, 0.0};
  const double a1 = M * q0 / q - 0.5 * q;
  const double a2 = M * q0 / q + 0.5 * q;
  // Re U_N = g M / (4 pi^2 q) [F(a1) - F(a2)], g = 4
  const double re = M / (kPi * kPi * q) * (LindhardF(a1, kF) - LindhardF(a2, kF));
  double im = 0;
  if (q0 > 0 && q0 < q) {
    const double w = std::min(2 * M * q0, kF * kF - a1 * a1);
    if (w > 0) im = -M / (2 * kPi * q) * w;
  }
  return {re, im};
}

std::complex<double> ValenciaRPA::LindhardN(double q0, double q, double kF, double M, double gamma) {
  if (gamma <= 0) return LindhardN(q0, q, kF, M);
  if (kF <= 0 || q <= 0) return {0.0, 0.0};
  // Direct + crossed ph bubbles: U_N(z) = g M / (4 pi^2 q) [G(a1) - G(a2)], a1,2 = M z / q -+ q / 2,
  // G(a) = a kF + (kF^2 - a^2)/2 ln[(a + kF)/(a - kF)], z = q0 + i gamma/2 (Im a > 0: the principal log
  // is continuous and reproduces the z -> q0 + i0 limit).
  const std::complex<double> z(q0, 0.5 * gamma);
  auto G = [kF](std::complex<double> a) { return a * kF + 0.5 * (kF * kF - a * a) * std::log((a + kF) / (a - kF)); };
  const std::complex<double> a1 = M * z / q - 0.5 * q, a2 = M * z / q + 0.5 * q;
  return M / (kPi * kPi * q) * (G(a1) - G(a2));
}

std::complex<double> ValenciaRPA::LindhardDelta(double q0, double q, double rho) const {
  const double omega = std::sqrt(par_.mDelta * par_.mDelta + q * q) - kMN;
  const double r = par_.fStarOverF;
  const std::complex<double> fwd = 1.0 / std::complex<double>(q0 - omega, 0.5 * par_.gammaDelta);
  return 4.0 / 9.0 * r * r * rho * (fwd - 1.0 / (q0 + omega));
}

double ValenciaRPA::Vl(double q0, double q) const {
  const double f2 = 4 * kPi * par_.f2Over4Pi;
  const double t = q0 * q0 - q * q;  // four-momentum squared (< 0 for QE)
  const double L2 = par_.lambdaPi * par_.lambdaPi, m2 = kMPi * kMPi;
  const double ff = (L2 - m2) / (L2 - t);
  return f2 / m2 * (ff * ff * q * q / (t - m2) + par_.gPrime);
}

double ValenciaRPA::Vt(double q0, double q) const {
  const double f2 = 4 * kPi * par_.f2Over4Pi;
  const double t = q0 * q0 - q * q;
  const double L2 = par_.lambdaRho * par_.lambdaRho, mr2 = par_.mRho * par_.mRho;
  const double ff = (L2 - mr2) / (L2 - t);
  return f2 / (kMPi * kMPi) * (par_.cRho * ff * ff * q * q / (t - mr2) + par_.gPrime);
}

ValenciaRPA::Factors ValenciaRPA::Evaluate(double q0In, double q, double rhoFm3, double rhoRatio) const {
  Factors f;
  if (rhoFm3 <= 0 || q <= 0 || q >= par_.qZero) return f;
  const double rho = rhoFm3 * kHbarC * kHbarC * kHbarC;  // GeV^3
  const double kF = std::cbrt(1.5 * kPi * kPi * rho);    // symmetric matter, rho = 2 kF^3 / (3 pi^2)
  // Clamp into the (open) 1p1h support of the local Fermi sea.
  const double lo = std::max(0.0, q * q - 2 * q * kF) / (2 * kMN), hi = (q * q + 2 * q * kF) / (2 * kMN);
  const double eps = 1e-3 * (hi - lo);
  const double q0 = std::clamp(q0In, lo + eps, hi - eps);
  const std::complex<double> UN = LindhardN(q0, q, kF, kMN, par_.gammaPH);
  const std::complex<double> U = UN + (par_.delta ? LindhardDelta(q0, q, rho) : std::complex<double>(0.0));
  const double x = std::clamp(rhoRatio, 0.0, 1.0);
  const double fp = x * par_.f0pIn + (1 - x) * par_.f0pEx;
  const double c0 = par_.c0 / (kHbarC * kHbarC * kHbarC);  // GeV^-2
  f.CN = 1.0 / std::norm(1.0 - c0 * fp * UN);
  f.CL = 1.0 / std::norm(1.0 - U * Vl(q0, q));
  f.CT = 1.0 / std::norm(1.0 - U * Vt(q0, q));
  if (q > par_.qFull) {
    const double s = (par_.qZero - q) / (par_.qZero - par_.qFull);
    f.CN = 1 + s * (f.CN - 1);
    f.CL = 1 + s * (f.CL - 1);
    f.CT = 1 + s * (f.CT - 1);
  }
  return f;
}

double ValenciaRPA::Contract(const FourVector& k, const FourVector& kp, const FourVector& p, const FourVector& pp,
                             const WeakFormFactors& ff, double M, int helicity, const Factors& fac) {
  const FourVector qv = pp - p;
  const FourVector Pv = pp + p;
  const double Q2 = -qv.M2();
  const double M2 = M * M;

  // Frame with z along the three-momentum transfer.
  const ThreeVector ez = qv.Vect().Unit();
  const ThreeVector a = std::abs(ez.x) < 0.9 ? ThreeVector{1, 0, 0} : ThreeVector{0, 1, 0};
  const ThreeVector ex = ez.Cross(a).Unit();
  const ThreeVector ey = ez.Cross(ex);
  double K[4], Kp[4], Pi[4], Pf[4], Pp[4], Qv[4];
  Components(k, ex, ey, ez, K);
  Components(kp, ex, ey, ez, Kp);
  Components(p, ex, ey, ez, Pi);
  Components(pp, ex, ey, ez, Pf);
  Components(Pv, ex, ey, ez, Pp);
  Components(qv, ex, ey, ez, Qv);
  const double g[4] = {1, -1, -1, -1};
  const double kk = k.Dot(kp);

  // Block-summed contractions  sum_{(m,n) in block} g_mm g_nn L^{mn} T^{mn}  of the symmetric lepton
  // tensor with the four symmetric hadronic structures S1, SA (p'p + pp' - g Y), PP and QQ.
  double bS1[kNBlocks] = {}, bSA[kNBlocks] = {}, bPP[kNBlocks] = {}, bQQ[kNBlocks] = {};
  for (int m = 0; m < 4; ++m) {
    for (int n = 0; n < 4; ++n) {
      const double gmn = m == n ? g[m] : 0.0;
      const double L = 8 * (K[m] * Kp[n] + Kp[m] * K[n] - gmn * kk);
      const double w = g[m] * g[n] * L;
      const Block b = BlockOf(m, n);
      const double pp2 = Pf[m] * Pi[n] + Pi[m] * Pf[n];
      bS1[b] += w * (pp2 - gmn * 0.5 * Q2);
      bSA[b] += w * (pp2 - gmn * (2 * M2 + 0.5 * Q2));
      bPP[b] += w * Pp[m] * Pp[n];
      bQQ[b] += w * Qv[m] * Qv[n];
    }
  }
  const double CN = fac.CN, CT = fac.CT, CL = fac.CL;
  // Factor per block for each form-factor structure (theory/local_fermi_gas.md, Table 2).
  const double f11[kNBlocks] = {CN, CN, CN, 1, 1, 1};
  const double f12[kNBlocks] = {CN, CN, CN, CT, 1, 1};
  const double f22[kNBlocks] = {1, 1, 1, CT, 1, 1};
  const double fAA[kNBlocks] = {1, CL, CL, CT, 1, 1};
  const double fP[kNBlocks] = {CL, CL, CL, 1, 1, 1};  // q^mu q^nu has only 00, 0z, zz components
  double c11 = 0, c12 = 0, c22 = 0, cAA = 0, cAP = 0, cPP = 0;
  const double kF2 = 0.5 * Q2 / M2 - 2;
  for (int b = 0; b < kNBlocks; ++b) {
    c11 += f11[b] * 4 * bS1[b];
    c12 += f12[b] * (8 * bS1[b] - 4 * bPP[b]);
    c22 += f22[b] * (4 * bS1[b] + kF2 * bPP[b]);
    cAA += fAA[b] * 4 * bSA[b];
    cAP += fP[b] * (-8) * bQQ[b];
    cPP += fP[b] * 2 * Q2 / M2 * bQQ[b];
  }
  const double F1 = ff.F1, F2 = ff.F2, FA = ff.FA, FP = ff.FP;
  const double sym = F1 * F1 * c11 + F1 * F2 * c12 + F2 * F2 * c22 + FA * FA * cAA + FA * FP * cAP + FP * FP * cPP;

  // Vector-axial interference: total (as in ContractQE) and its transverse (xy) part, which gets C_T.
  const double vaCoef = helicity * (F1 + F2) * FA;
  const double anti = 128 * vaCoef * (k.Dot(p) * kp.Dot(qv) - k.Dot(qv) * kp.Dot(p));
  const double antiTT = -128 * vaCoef * (Qv[0] * Pi[3] - Qv[3] * Pi[0]) * (Kp[0] * K[3] - Kp[3] * K[0]);
  return sym + anti + (CT - 1) * antiTT;
}

}  // namespace nuflow
