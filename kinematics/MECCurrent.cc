#include "kinematics/MECCurrent.h"

#include <cmath>
#include <sstream>

#include "core/Constants.h"

namespace nuflow {

using namespace constants;
using dirac::Comp;
using dirac::Gammas;
using dirac::Mat;
using dirac::Metric;
using cplx = std::complex<double>;

namespace {

constexpr double kM = kMassNucleon;
constexpr double kMpi = kMassPionAverage;

// Dirac spinor u(p, s) (Dirac representation, ubar u = 2M), spin s = 0 (up) / 1 (down) along z, and its
// Dirac adjoint as a row vector.
struct Spinor {
  cplx c[4];
};

Spinor U(const FourVector& p, int s) {
  const double n = std::sqrt(p.e + kM);
  const double inv = 1.0 / (p.e + kM);
  const cplx I(0, 1);
  Spinor u;
  if (s == 0) {
    u.c[0] = n;
    u.c[1] = 0;
    u.c[2] = n * p.pz * inv;
    u.c[3] = n * (p.px + I * p.py) * inv;
  } else {
    u.c[0] = 0;
    u.c[1] = n;
    u.c[2] = n * (p.px - I * p.py) * inv;
    u.c[3] = -n * p.pz * inv;
  }
  return u;
}

Spinor Bar(const Spinor& u) { return {{std::conj(u.c[0]), std::conj(u.c[1]), -std::conj(u.c[2]), -std::conj(u.c[3])}}; }

// out[s'][s] = ubar(s') O u(s)
void Sandwich(const Mat& O, const Spinor ub[2], const Spinor u[2], cplx out[2][2]) {
  for (int s = 0; s < 2; ++s) {
    cplx Ou[4];
    for (int a = 0; a < 4; ++a) {
      cplx t = 0;
      for (int b = 0; b < 4; ++b) t += O.m[4 * a + b] * u[s].c[b];
      Ou[a] = t;
    }
    for (int sp = 0; sp < 2; ++sp) {
      cplx t = 0;
      for (int a = 0; a < 4; ++a) t += ub[sp].c[a] * Ou[a];
      out[sp][s] = t;
    }
  }
}

// One nucleon line i -> f: all bilinears needed by the MEC.
struct Line {
  Spinor ub[2], u[2];
  cplx g5[2][2];           // ubar gamma5 u
  cplx gam[4][2][2];       // ubar gamma^mu u
  cplx g5gam[4][2][2];     // ubar gamma5 gamma^mu u
  cplx qs[2][2];           // ubar Qslash u
  cplx dF[4][2][2];        // ubar O_forward^mu u (Delta excited on this line)
  cplx dB[4][2][2];        // ubar O_backward^mu u
};

double Monopole(double lambda, double k2) {
  if (lambda <= 0) return 1.0;
  const double L2 = lambda * lambda;
  return (L2 - kMpi * kMpi) / (L2 - k2);
}

}  // namespace

std::string MECCurrent::Name() const {
  std::ostringstream s;
  s.precision(6);
  s << "MEC-RS17[C5A0=" << par_.C5A0 << ",MAD=" << par_.MADelta << ",C3V0=" << par_.C3V0 << ",f*=" << par_.fStar
    << ",ReD=" << (par_.realDeltaPropagator ? 1 : 0) << ",LpiNN=" << par_.lambdaPiNN << ",LpiND=" << par_.lambdaPiNDelta
    << ",A=" << par_.axialScale << ",V=" << par_.vectorScale << "," << nucleonFF_.Name() << "]";
  return s.str();
}

double MECCurrent::DeltaWidth(double s) const {
  if (s <= (kM + kMpi) * (kM + kMpi)) return 0;
  const double rs = std::sqrt(s);
  const double lam = (s - (kM + kMpi) * (kM + kMpi)) * (s - (kM - kMpi) * (kM - kMpi));
  const double pcm = std::sqrt(std::max(0.0, lam)) / (2 * rs);
  const double f = par_.fStar / kMpi;
  return f * f / (6 * kPi) * (kM / rs) * pcm * pcm * pcm;
}

void MECCurrent::ComputeCurrents(const FourVector& Q, const FourVector& f1, const FourVector& i1, const FourVector& f2,
                                 const FourVector& i2, Currents& out) const {
  const Gammas& G = Gammas::Get();
  const double MD = par_.MDelta, MD2 = MD * MD;
  const double Q2 = -Q.M2();  // > 0 for spacelike Q
  const double mpi2 = kMpi * kMpi;

  // Form factors.
  const WeakFormFactors nff = nucleonFF_.Evaluate(Q2);
  const double F1V = nff.F1 * par_.vectorScale;
  const double MV2 = par_.MV * par_.MV;
  const double dV = 1 + Q2 / MV2;
  const double C3V = par_.vectorScale * par_.C3V0 / (dV * dV) / (1 + Q2 / (4 * MV2));
  const double dA = 1 + Q2 / (par_.MADelta * par_.MADelta);
  const double C5A = par_.axialScale * par_.C5A0 / (dA * dA);
  const double fm = par_.gA / (2 * par_.fPi);  // f_piNN / m_pi
  const double fsm = par_.fStar / kMpi;        // f* / m_pi
  const double ax = par_.axialScale / par_.gA;

  const FourVector k11 = f1 - i1, k22 = f2 - i2;
  const double t11 = k11.M2(), t22 = k22.M2();
  const double Frho11 = 1 / (1 - t11 / (par_.mRho * par_.mRho)), Frho22 = 1 / (1 - t22 / (par_.mRho * par_.mRho));
  const double FNN11 = Monopole(par_.lambdaPiNN, t11), FNN22 = Monopole(par_.lambdaPiNN, t22);
  const double FND11 = Monopole(par_.lambdaPiNDelta, t11), FND22 = Monopole(par_.lambdaPiNDelta, t22);

  const Mat Qs = G.Slash(Q);
  const Mat Qs5 = Qs * G.g5;

  // Delta operators O^mu for the forward (P = i + Q) and backward (P = f - Q) terms, contracted with the
  // pion momentum k that is absorbed by the other line (RS17 Eqs. 33-40, C3V and C5A only).
  auto deltaProp = [&](double s) -> cplx {
    const double d = s - MD2, g = MD * DeltaWidth(s);
    if (par_.realDeltaPropagator) return d / (d * d + g * g);
    return 1.0 / cplx(d, g);
  };
  auto forward = [&](const FourVector& P, const FourVector& k, Mat O[4]) {
    const cplx D = deltaProp(P.M2());
    const Mat pre = (G.Slash(P) + G.one * MD) * (-1.0);
    const Mat ks = G.Slash(k);
    const double kP = k.Dot(P), QP = Q.Dot(P), kQ = k.Dot(Q);
    Mat B = G.one * (kQ - 2 * kP * QP / (3 * MD2));
    B.AddScaled(ks * Qs, -1.0 / 3);
    B.AddScaled(Qs, kP / (3 * MD));
    B.AddScaled(ks, -QP / (3 * MD));
    const Mat Xq = pre * B;
    Mat R = Qs5 * (C3V / kM);
    R.AddScaled(G.one, C5A);
    for (int mu = 0; mu < 4; ++mu) {
      Mat Bm = G.one * (Comp(k, mu) - 2 * kP * Comp(P, mu) / (3 * MD2));
      Bm.AddScaled(ks * G.g[mu], -1.0 / 3);
      Bm.AddScaled(G.g[mu], kP / (3 * MD));
      Bm.AddScaled(ks, -Comp(P, mu) / (3 * MD));
      const Mat X = pre * Bm;
      O[mu] = (X * R + Xq * (G.g5g[mu] * (-C3V / kM))) * D;
    }
  };
  auto backward = [&](const FourVector& P, const FourVector& k, Mat O[4]) {
    const cplx D = deltaProp(P.M2());
    const Mat pre = (G.Slash(P) + G.one * MD) * (-1.0);
    const Mat ks = G.Slash(k);
    const double kP = k.Dot(P), QP = Q.Dot(P), kQ = k.Dot(Q);
    Mat B = G.one * (kQ - 2 * QP * kP / (3 * MD2));
    B.AddScaled(Qs * ks, -1.0 / 3);
    B.AddScaled(ks, QP / (3 * MD));
    B.AddScaled(Qs, -kP / (3 * MD));
    const Mat qY = pre * B;
    Mat L = Qs5 * (-C3V / kM);
    L.AddScaled(G.one, C5A);
    for (int mu = 0; mu < 4; ++mu) {
      Mat Bm = G.one * (Comp(k, mu) - 2 * Comp(P, mu) * kP / (3 * MD2));
      Bm.AddScaled(G.g[mu] * ks, -1.0 / 3);
      Bm.AddScaled(ks, Comp(P, mu) / (3 * MD));
      Bm.AddScaled(G.g[mu], -kP / (3 * MD));
      const Mat Y = pre * Bm;
      O[mu] = (L * Y + (G.g5g[mu] * (C3V / kM)) * qY) * D;
    }
  };

  Line L1, L2;
  auto fill = [&](Line& L, const FourVector& f, const FourVector& i, const FourVector& kOther) {
    for (int s = 0; s < 2; ++s) {
      L.u[s] = U(i, s);
      L.ub[s] = Bar(U(f, s));
    }
    Sandwich(G.g5, L.ub, L.u, L.g5);
    Sandwich(Qs, L.ub, L.u, L.qs);
    for (int mu = 0; mu < 4; ++mu) {
      Sandwich(G.g[mu], L.ub, L.u, L.gam[mu]);
      Sandwich(G.g5 * G.g[mu], L.ub, L.u, L.g5gam[mu]);
    }
    Mat O[4];
    forward(i + Q, kOther, O);
    for (int mu = 0; mu < 4; ++mu) Sandwich(O[mu], L.ub, L.u, L.dF[mu]);
    backward(f - Q, kOther, O);
    for (int mu = 0; mu < 4; ++mu) Sandwich(O[mu], L.ub, L.u, L.dB[mu]);
  };
  fill(L1, f1, i1, k22);  // Delta on line 1 emits the pion k22 absorbed by line 2
  fill(L2, f2, i2, k11);

  // V_piNN(p', h) = ubar gamma5 kslash u / (k^2 - m_pi^2) = -2M ubar gamma5 u / (k^2 - m_pi^2)  (RS17 Eqs. 21-22)
  const double prop11 = -2 * kM / (t11 - mpi2) * FNN11, prop22 = -2 * kM / (t22 - mpi2) * FNN22;
  const double cDelta = -fsm * fm;
  const double r6 = 1.0 / std::sqrt(6.0);
  const double Qpole = 1.0 / (Q.M2() - mpi2);

  for (int s = 0; s < 16; ++s) {
    const int sf1 = s & 1, sf2 = (s >> 1) & 1, si1 = (s >> 2) & 1, si2 = (s >> 3) & 1;
    const cplx V1 = prop11 * L1.g5[sf1][si1];
    const cplx V2 = prop22 * L2.g5[sf2][si2];
    for (int mu = 0; mu < 4; ++mu) {
      // Seagull (RS17 Eqs. 17-18); the contact vertex carries a second piNN form factor when enabled.
      const cplx sea = fm * fm *
                       (F1V * (V1 * FNN11 * L2.g5gam[mu][sf2][si2] - V2 * FNN22 * L1.g5gam[mu][sf1][si1]) +
                        ax * (Frho22 * V1 * FNN11 * L2.gam[mu][sf2][si2] - Frho11 * V2 * FNN22 * L1.gam[mu][sf1][si1]));
      // Pion in flight (RS17 Eq. 26), vector only.
      const cplx pif = fm * fm * F1V * V1 * V2 * (Comp(k11, mu) - Comp(k22, mu));
      // Pion pole (RS17 Eq. 31), axial only.
      const cplx pole = fm * fm * ax * Comp(Q, mu) * Qpole *
                        (Frho11 * V2 * L1.qs[sf1][si1] - Frho22 * V1 * L2.qs[sf2][si2]);
      // Delta forward/backward with the Delta on line 1 (F1, B1) or line 2 (F2, B2) (RS17 Eqs. 33-34).
      const cplx F1 = cDelta * FND22 * V2 * L1.dF[mu][sf1][si1];
      const cplx B1 = cDelta * FND22 * V2 * L1.dB[mu][sf1][si1];
      const cplx F2 = cDelta * FND11 * V1 * L2.dF[mu][sf2][si2];
      const cplx B2 = cDelta * FND11 * V1 * L2.dB[mu][sf2][si2];
      // Isospin decomposition, RS17 Eqs. (45)-(53).
      out.J1[s][mu] = -2 * r6 * (F2 + B2);
      out.J2[s][mu] = -2 * r6 * (F1 + B1);
      out.J3[s][mu] = par_.backgroundScale * (sea + pif + pole) + r6 * (F1 - B1 - F2 + B2);
    }
  }
}

void MECCurrent::Antisymmetrized(MECChannel ch, const FourVector& Q, const FourVector& h1, const FourVector& h2,
                                 const FourVector& p1, const FourVector& p2, Amplitudes& out) const {
  Currents d, x;
  ComputeCurrents(Q, p1, h1, p2, h2, d);
  const bool finalExchange = (ch == MECChannel::kNuNPtoPP || ch == MECChannel::kNubarNPtoNN);
  // Exchange term: final momenta swapped (p2' on line 1) or initial momenta swapped (h2 on line 1).
  if (finalExchange) ComputeCurrents(Q, p2, h1, p1, h2, x);
  else ComputeCurrents(Q, p1, h2, p2, h1, x);
  auto eff = [&](const Currents& c, int s, int mu) -> cplx {
    switch (ch) {
      case MECChannel::kNuNPtoPP: return c.J1[s][mu] + c.J3[s][mu];
      case MECChannel::kNuNNtoNP: return c.J2[s][mu] + c.J3[s][mu];
      case MECChannel::kNubarPPtoNP: return c.J1[s][mu] - c.J3[s][mu];
      case MECChannel::kNubarNPtoNN: return c.J2[s][mu] - c.J3[s][mu];
    }
    return 0;
  };
  for (int s = 0; s < 16; ++s) {
    const int sf1 = s & 1, sf2 = (s >> 1) & 1, si1 = (s >> 2) & 1, si2 = (s >> 3) & 1;
    // Spin index of the exchanged labelling.
    const int sx = finalExchange ? (sf2 | (sf1 << 1) | (si1 << 2) | (si2 << 3))
                                 : (sf1 | (sf2 << 1) | (si2 << 2) | (si1 << 3));
    for (int mu = 0; mu < 4; ++mu) out.a[s][mu] = eff(d, s, mu) - eff(x, sx, mu);
  }
}

void LeptonTensor(const FourVector& k, const FourVector& kp, double mLep, int helicity, cplx L[4][4]) {
  const Gammas& G = Gammas::Get();
  const Mat PL = G.one - G.g5;
  const Mat Kin = helicity > 0 ? G.Slash(k) : G.Slash(kp) - G.one * mLep;
  const Mat Kout = helicity > 0 ? G.Slash(kp) + G.one * mLep : G.Slash(k);
  Mat gPL[4], Um[4];
  for (int mu = 0; mu < 4; ++mu) gPL[mu] = G.g[mu] * PL;
  for (int mu = 0; mu < 4; ++mu) Um[mu] = Kout * gPL[mu] * Kin;
  for (int mu = 0; mu < 4; ++mu)
    for (int nu = 0; nu < 4; ++nu) L[mu][nu] = dirac::TraceProduct(Um[mu], gPL[nu]);
}

double ContractLepton(const cplx L[4][4], const cplx A[4]) {
  cplx Al[4];
  for (int mu = 0; mu < 4; ++mu) Al[mu] = Metric(mu) * A[mu];
  cplx t = 0;
  for (int mu = 0; mu < 4; ++mu)
    for (int nu = 0; nu < 4; ++nu) t += L[mu][nu] * Al[mu] * std::conj(Al[nu]);
  return t.real();
}

void MECHadronTensor(const MECCurrent& cur, MECChannel ch, double q, double omegaEff, double kF, double Z, long n,
                     Random& rng, cplx W[4][4], bool directOnly) {
  for (int mu = 0; mu < 4; ++mu)
    for (int nu = 0; nu < 4; ++nu) W[mu][nu] = 0;
  const FourVector Q(omegaEff, 0, 0, q);
  auto inSphere = [&]() {
    const double k = kF * std::cbrt(rng.UniformOpen());
    const double c = 2 * rng.Uniform() - 1, s = std::sqrt(std::max(0.0, 1 - c * c)), ph = 2 * kPi * rng.Uniform();
    const ThreeVector v{k * s * std::cos(ph), k * s * std::sin(ph), k * c};
    return FourVector(std::sqrt(v.Mag2() + kM * kM), v);
  };
  for (long i = 0; i < n; ++i) {
    const FourVector h1 = inSphere(), h2 = inSphere();
    const FourVector P = Q + h1 + h2;
    const double s = P.M2();
    if (P.e <= 0 || s <= 4 * kM * kM) continue;
    const double rs = std::sqrt(s), pstar = std::sqrt(s / 4 - kM * kM);
    const double c = 2 * rng.Uniform() - 1, sn = std::sqrt(std::max(0.0, 1 - c * c)), ph = 2 * kPi * rng.Uniform();
    const FourVector p1 = FourVector(rs / 2, pstar * sn * std::cos(ph), pstar * sn * std::sin(ph), pstar * c).Boosted(P.BoostVector());
    const FourVector p2 = P - p1;
    if (p1.P2() < kF * kF || p2.P2() < kF * kF) continue;
    // int d3p1' delta(E)/(E1'E2') = int dOmega* p*/W;  M^4/(E...) r(ubar u = 1) = r(covariant)/(16 E1 E2 E1' E2')
    const double jac = 4 * kPi * pstar / rs / (16 * h1.e * h2.e);
    MECCurrent::Amplitudes A;
    double fac = 2;
    if (directOnly) {
      MECCurrent::Currents d;
      cur.ComputeCurrents(Q, p1, h1, p2, h2, d);
      for (int sp = 0; sp < 16; ++sp)
        for (int mu = 0; mu < 4; ++mu) {
          switch (ch) {
            case MECChannel::kNuNPtoPP: A.a[sp][mu] = d.J1[sp][mu] + d.J3[sp][mu]; break;
            case MECChannel::kNuNNtoNP: A.a[sp][mu] = d.J2[sp][mu] + d.J3[sp][mu]; break;
            case MECChannel::kNubarPPtoNP: A.a[sp][mu] = d.J1[sp][mu] - d.J3[sp][mu]; break;
            case MECChannel::kNubarNPtoNN: A.a[sp][mu] = d.J2[sp][mu] - d.J3[sp][mu]; break;
          }
        }
      fac = 4;
    } else {
      cur.Antisymmetrized(ch, Q, h1, h2, p1, p2, A);
    }
    for (int sp = 0; sp < 16; ++sp)
      for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu) W[mu][nu] += fac * jac * std::conj(A.a[sp][mu]) * A.a[sp][nu];
  }
  const double V = 3 * kPi * kPi * Z / (kF * kF * kF);
  const double sphere = 4 * kPi / 3 * kF * kF * kF;
  const double pref = V * sphere * sphere / std::pow(2 * kPi, 9) / static_cast<double>(n);
  for (int mu = 0; mu < 4; ++mu)
    for (int nu = 0; nu < 4; ++nu) W[mu][nu] *= pref;
}

}  // namespace nuflow
