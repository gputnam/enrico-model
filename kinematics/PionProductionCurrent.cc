#include "kinematics/PionProductionCurrent.h"

#include <cmath>
#include <sstream>

#include "core/Constants.h"

namespace nuflow {

using namespace constants;
using dirac::cplx;
using dirac::Comp;
using dirac::Gammas;
using dirac::Mat;
using dirac::Metric;

std::string PionProductionCurrent::Name() const {
  std::ostringstream s;
  s.precision(6);
  s << "HNV[C5A0=" << par_.C5A0 << ",MAD=" << par_.MADelta << ",bg=" << (par_.background ? 1 : 0) << ","
    << nucleonFF_.Name() << "]";
  return s.str();
}

double PionProductionCurrent::DeltaWidth(double s) const {
  const double M = kMassNucleon, mpi = kMassPionAverage;
  if (s <= (M + mpi) * (M + mpi)) return 0;
  const double rs = std::sqrt(s);
  const double lam = (s - (M + mpi) * (M + mpi)) * (s - (M - mpi) * (M - mpi));
  const double pcm = std::sqrt(std::max(0.0, lam)) / (2 * rs);
  const double f = par_.fStar / mpi;
  return f * f / (6 * kPi) * (M / rs) * pcm * pcm * pcm;
}

NDeltaFormFactors PionProductionCurrent::DeltaFormFactors(double Q2) const {
  // HNV Eq. (46) (vector, Lalakulich-Paschos-Piranishvili), dipole C5A (arXiv:1001.4416) with Adler's
  // constraints C3A = 0, C4A = -C5A/4 and PCAC C6A = C5A M^2 / (mpi^2 + Q2) (HNV Eqs. 47-48).
  const double M2 = kMassNucleon * kMassNucleon, mpi2 = kMassPionAverage * kMassPionAverage;
  const double MV2 = par_.MV * par_.MV;
  const double dV = 1 + Q2 / MV2, gV = 1 / (dV * dV);
  NDeltaFormFactors c;
  c.C3V = par_.C3V0 * gV / (1 + Q2 / (4 * MV2));
  c.C4V = par_.C4V0 * gV / (1 + Q2 / (4 * MV2));
  c.C5V = par_.C5V0 * gV / (1 + Q2 / (0.776 * MV2));
  const double dA = 1 + Q2 / (par_.MADelta * par_.MADelta);
  c.C5A = par_.axialScale * par_.C5A0 / (dA * dA);
  c.C4A = -c.C5A / 4;
  c.C3A = 0;
  c.C6A = c.C5A * M2 / (mpi2 + Q2);
  return c;
}

void PionProductionCurrent::Pieces(const FourVector& p, const FourVector& pp, const FourVector& kpi,
                                   PionCurrentPieces& out) const {
  const Gammas& G = Gammas::Get();
  const cplx I(0, 1);
  const double M = kMassNucleon, M2 = M * M, mpi = kMassPionAverage, mpi2 = mpi * mpi;
  const double MD = par_.MDelta, MD2 = MD * MD;
  const FourVector q = pp + kpi - p;
  const double q2 = q.M2();
  const double Q2 = -q2;

  const NDeltaFormFactors C = DeltaFormFactors(Q2);
  const double C3V = C.C3V, C4V = C.C4V, C5V = C.C5V, C3A = C.C3A, C4A = C.C4A, C5A = C.C5A, C6A = C.C6A;

  const Mat qs = G.Slash(q), ks = G.Slash(kpi);
  const Mat qs5 = qs * G.g5, ks5 = ks * G.g5;
  const double kq = kpi.Dot(q);

  // ---- Direct Delta pole: i (f*/mpi) sqrt3 kpi^a P_ab(pD) Gamma^{b mu}(p,q) / (pD^2 - MD^2 + i MD Gamma)
  {
    const FourVector pD = p + q;
    const double s = pD.M2();
    const cplx fac = I * (par_.fStar / mpi) * std::sqrt(3.0) / cplx(s - MD2, MD * DeltaWidth(s));
    const Mat pre = (G.Slash(pD) + G.one * MD) * (-1.0);
    const double kpD = kpi.Dot(pD), qpD = q.Dot(pD), qp = q.Dot(p);
    // X^mu = kpi_a P^{a mu}, Xq = kpi_a P^{a b} q_b
    Mat Xq;
    {
      Mat B = G.one * (kq - 2 * kpD * qpD / (3 * MD2));
      B.AddScaled(ks * qs, -1.0 / 3);
      B.AddScaled(qs, kpD / (3 * MD));
      B.AddScaled(ks, -qpD / (3 * MD));
      Xq = pre * B;
    }
    const double a1 = (C4V * qpD + C5V * qp) / M2;
    const double a2 = C4A * qpD / M2 + C5A;
    Mat R = qs5 * (C3V / M);
    R.AddScaled(G.g5, a1);
    R.AddScaled(G.one, a2);
    R.AddScaled(qs, C3A / M);
    for (int mu = 0; mu < 4; ++mu) {
      Mat B = G.one * (Comp(kpi, mu) - 2 * kpD * Comp(pD, mu) / (3 * MD2));
      B.AddScaled(ks * G.g[mu], -1.0 / 3);
      B.AddScaled(G.g[mu], kpD / (3 * MD));
      B.AddScaled(ks, -Comp(pD, mu) / (3 * MD));
      const Mat X = pre * B;
      const double b1 = (C4V * Comp(pD, mu) + C5V * Comp(p, mu)) / M2;
      const double b2 = (-C4A * Comp(pD, mu) + C6A * Comp(q, mu)) / M2;
      Mat T = G.g5g[mu] * (-C3V / M);
      T.AddScaled(G.g5, -b1);
      T.AddScaled(G.one, b2);
      T.AddScaled(G.g[mu], -C3A / M);
      out.dp[mu] = (X * R + Xq * T) * fac;
    }
  }

  if (!par_.background) {
    for (int mu = 0; mu < 4; ++mu) out.cdp[mu] = out.np[mu] = out.cnp[mu] = out.ct[mu] = out.pp[mu] = out.pf[mu] = Mat{};
    return;
  }

  // ---- Crossed Delta pole: i (f*/mpi)/sqrt3 Gammahat^{mu a}(p',q) P_ab(pD) kpi^b / (pD^2 - MD^2 + i MD Gamma),
  //      pD = p' - q, Gammahat^{mu a}(p',q) = gamma0 [Gamma^{a mu}(p',-q)]^dagger gamma0.
  {
    const FourVector pD = pp - q;
    const double u = pD.M2();
    const cplx fac = I * (par_.fStar / mpi) / std::sqrt(3.0) / cplx(u - MD2, MD * DeltaWidth(u));
    const Mat pre = (G.Slash(pD) + G.one * MD) * (-1.0);
    const double kpD = kpi.Dot(pD), qpD = q.Dot(pD), qpp = q.Dot(pp);
    Mat qY;  // q^a P_ab kpi^b
    {
      Mat B = G.one * (kq - 2 * qpD * kpD / (3 * MD2));
      B.AddScaled(qs * ks, -1.0 / 3);
      B.AddScaled(ks, qpD / (3 * MD));
      B.AddScaled(qs, -kpD / (3 * MD));
      qY = pre * B;
    }
    const double a1 = (C4V * qpD + C5V * qpp) / M2;
    Mat L = qs5 * (-C3V / M);  // multiplies Y^mu from the left
    L.AddScaled(G.g5, a1);
    L.AddScaled(G.one, C5A - C4A * qpD / M2);
    L.AddScaled(qs, -C3A / M);
    for (int mu = 0; mu < 4; ++mu) {
      Mat B = G.one * (Comp(kpi, mu) - 2 * Comp(pD, mu) * kpD / (3 * MD2));
      B.AddScaled(G.g[mu] * ks, -1.0 / 3);
      B.AddScaled(ks, Comp(pD, mu) / (3 * MD));
      B.AddScaled(G.g[mu], -kpD / (3 * MD));
      const Mat Y = pre * B;
      const double b1 = (C4V * Comp(pD, mu) + C5V * Comp(pp, mu)) / M2;
      Mat T = G.g5g[mu] * (C3V / M);
      T.AddScaled(G.g5, -b1);
      T.AddScaled(G.one, (C4A * Comp(pD, mu) + C6A * Comp(q, mu)) / M2);
      T.AddScaled(G.g[mu], C3A / M);
      out.cdp[mu] = (L * Y + T * qY) * fac;
    }
  }

  // ---- Chiral background (HNV Eq. 51), nucleon form factors from NucleonFormFactors (F1, F2 isovector;
  //      FA plays the role of G_A).
  WeakFormFactors ff = nucleonFF_.Evaluate(Q2);
  ff.FA *= par_.axialScale;
  const double cBg = par_.gA / (std::sqrt(2.0) * par_.fPi);
  const double t = (q - kpi).M2();
  const double Frho = 1 / (1 - t / (par_.mRho * par_.mRho));
  const double FrhoA = par_.axialScale * Frho;  // axial pieces: CT (F_rho part) and pion pole
  const FourVector pS = p + q, pU = pp - q;
  const double sN = pS.M2() - M2, uN = pU.M2() - M2;
  const Mat preNP = ks5 * (G.Slash(pS) + G.one * M) * (1.0 / sN);
  const Mat postCNP = (G.Slash(pU) + G.one * M) * ks5 * (1.0 / uN);
  const Mat ctCore = G.g5 * (par_.gA * ff.F1) - G.one * FrhoA;
  for (int mu = 0; mu < 4; ++mu) {
    // V_N^mu - A_N^mu  (HNV Eq. 34), sigma^{mu nu} q_nu = (i/2)[gamma^mu, qslash]
    Mat VA = G.g[mu] * ff.F1;
    VA.AddScaled(G.g[mu] * qs - qs * G.g[mu], -ff.F2 / (4 * M));
    VA.AddScaled(G.g5g[mu], -ff.FA);
    VA.AddScaled(qs5, -ff.FA * Comp(q, mu) / (mpi2 - q2));
    out.np[mu] = preNP * VA * (-I * cBg);
    out.cnp[mu] = VA * postCNP * (-I * cBg);
    out.ct[mu] = G.g[mu] * ctCore * (-I / (std::sqrt(2.0) * par_.fPi));
    out.pp[mu] = qs * (-I * FrhoA / (std::sqrt(2.0) * par_.fPi) * Comp(q, mu) / (q2 - mpi2));
    out.pf[mu] = G.g5 * (-I * (ff.F1 * cBg * (2 * Comp(kpi, mu) - Comp(q, mu)) / (t - mpi2) * 2 * M));
  }
}

void PionProductionCurrent::Combine(const PionCurrentPieces& P, PionIsospin iso, Mat S[4]) const {
  // Coefficients (C_DP, C_CDP, C_NP, C_CNP, C_CT, C_PP, C_PF) of HNV Eq. (51).
  static const double r2 = std::sqrt(0.5);
  double c[7] = {0, 0, 0, 0, 0, 0, 0};
  switch (iso) {
    case PionIsospin::PPiPlus: { const double a[7] = {1, 1, 0, 1, 1, 1, 1}; std::copy(a, a + 7, c); break; }
    case PionIsospin::NPiPlus: { const double a[7] = {1.0 / 3, 3, 1, 0, -1, -1, -1}; std::copy(a, a + 7, c); break; }
    case PionIsospin::PPiZero: {
      // -(PPiPlus - NPiPlus)/sqrt2
      const double a[7] = {-(2.0 / 3) * r2, 2 * r2, r2, -r2, -2 * r2, -2 * r2, -2 * r2};
      std::copy(a, a + 7, c);
      break;
    }
  }
  for (int mu = 0; mu < 4; ++mu) {
    Mat s = P.dp[mu] * c[0];
    s.AddScaled(P.cdp[mu], c[1]);
    s.AddScaled(P.np[mu], c[2]);
    s.AddScaled(P.cnp[mu], c[3]);
    s.AddScaled(P.ct[mu], c[4]);
    s.AddScaled(P.pp[mu], c[5]);
    s.AddScaled(P.pf[mu], c[6]);
    S[mu] = s;
  }
}

double PionProductionCurrent::Contract(const FourVector& k, const FourVector& kp, double mLep, int helicity,
                                       const FourVector& p, double mIn, const FourVector& pp, double mOut,
                                       const Mat S[4]) {
  const Gammas& G = Gammas::Get();
  static const Mat PL = Gammas::Get().one - Gammas::Get().g5;
  static const Mat gPL[4] = {G.g[0] * PL, G.g[1] * PL, G.g[2] * PL, G.g[3] * PL};
  // Lepton line: nu: Tr[(k'slash + m) g^mu PL kslash g^nu PL]; nubar: Tr[kslash g^mu PL (k'slash - m) g^nu PL].
  const Mat Kin = helicity > 0 ? G.Slash(k) : G.Slash(kp) - G.one * mLep;
  const Mat Kout = helicity > 0 ? G.Slash(kp) + G.one * mLep : G.Slash(k);
  const Mat Pp = G.Slash(pp) + G.one * mOut, Pi = G.Slash(p) + G.one * mIn;
  Mat U[4], A[4], B[4];
  for (int mu = 0; mu < 4; ++mu) {
    U[mu] = Kout * gPL[mu] * Kin;
    A[mu] = Pp * S[mu] * Pi;
    B[mu] = G.Bar(S[mu]);
  }
  double total = 0;
  for (int mu = 0; mu < 4; ++mu)
    for (int nu = 0; nu < 4; ++nu) {
      const cplx L = dirac::TraceProduct(U[mu], gPL[nu]);
      const cplx H = dirac::TraceProduct(A[mu], B[nu]);
      total += (L * H).real() * Metric(mu) * Metric(nu);
    }
  return total;
}

}  // namespace nuflow
