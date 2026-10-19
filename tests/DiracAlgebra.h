#pragma once
// Brute-force Dirac algebra (explicit 4x4 gamma matrices in the Dirac representation), used by the
// tests to validate the closed-form QE contraction independently of any trace identities.

#include <array>
#include <complex>

#include "core/FourVector.h"
#include "kinematics/FormFactors.h"

namespace nuflow::test {

using cplx = std::complex<double>;
struct Mat {
  std::array<cplx, 16> m{};
  cplx& operator()(int i, int j) { return m[4 * i + j]; }
  cplx operator()(int i, int j) const { return m[4 * i + j]; }
  Mat operator*(const Mat& o) const {
    Mat r;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        cplx s = 0;
        for (int k = 0; k < 4; ++k) s += (*this)(i, k) * o(k, j);
        r(i, j) = s;
      }
    return r;
  }
  Mat operator+(const Mat& o) const {
    Mat r;
    for (int i = 0; i < 16; ++i) r.m[i] = m[i] + o.m[i];
    return r;
  }
  Mat operator-(const Mat& o) const {
    Mat r;
    for (int i = 0; i < 16; ++i) r.m[i] = m[i] - o.m[i];
    return r;
  }
  Mat operator*(cplx s) const {
    Mat r;
    for (int i = 0; i < 16; ++i) r.m[i] = m[i] * s;
    return r;
  }
  Mat Dagger() const {
    Mat r;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) r(i, j) = std::conj((*this)(j, i));
    return r;
  }
  cplx Trace() const { return m[0] + m[5] + m[10] + m[15]; }
  static Mat Identity() {
    Mat r;
    for (int i = 0; i < 4; ++i) r(i, i) = 1;
    return r;
  }
};

struct Gammas {
  Mat g[4], g5, one;
  Gammas() {
    const cplx I(0, 1);
    one = Mat::Identity();
    g[0](0, 0) = 1; g[0](1, 1) = 1; g[0](2, 2) = -1; g[0](3, 3) = -1;
    // gamma^i = [[0, sigma_i], [-sigma_i, 0]]
    const cplx s[3][2][2] = {{{0, 1}, {1, 0}}, {{0, -I}, {I, 0}}, {{1, 0}, {0, -1}}};
    for (int i = 0; i < 3; ++i)
      for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) {
          g[i + 1](a, b + 2) = s[i][a][b];
          g[i + 1](a + 2, b) = -s[i][a][b];
        }
    g5 = g[0] * g[1] * g[2] * g[3] * I;
  }
  Mat Slash(const FourVector& v) const { return g[0] * v.e - g[1] * v.px - g[2] * v.py - g[3] * v.pz; }
  Mat Bar(const Mat& x) const { return g[0] * x.Dagger() * g[0]; }
};

inline double Metric(int mu) { return mu == 0 ? 1.0 : -1.0; }
inline double Comp(const FourVector& v, int mu) {
  return mu == 0 ? v.e : mu == 1 ? v.px : mu == 2 ? v.py : v.pz;
}

// L_{mu nu} H^{mu nu} computed by explicit traces. axialSign multiplies the whole axial current
// (F_A and F_P terms); the Standard-Model V-A nucleon current corresponds to axialSign = -1 for F_A > 0.
inline double BruteForceContraction(const FourVector& k, const FourVector& kp, double mLep, const FourVector& p,
                                    const FourVector& pp, double M, const WeakFormFactors& ff, int helicity,
                                    int axialSign) {
  static const Gammas G;
  const cplx I(0, 1);
  const FourVector q = pp - p;
  const Mat PL = G.one - G.g5;
  Mat Gam[4], GamBar[4];
  for (int mu = 0; mu < 4; ++mu) {
    Mat sigq;  // sigma^{mu alpha} q_alpha
    for (int al = 0; al < 4; ++al) {
      const Mat sig = (G.g[mu] * G.g[al] - G.g[al] * G.g[mu]) * (0.5 * I);
      sigq = sigq + sig * (Metric(al) * Comp(q, al));
    }
    Gam[mu] = G.g[mu] * ff.F1 + sigq * (I * ff.F2 / (2 * M)) + G.g[mu] * G.g5 * (axialSign * ff.FA) +
              G.g5 * (axialSign * Comp(q, mu) * ff.FP / M);
    GamBar[mu] = G.Bar(Gam[mu]);
  }
  const Mat Pp = G.Slash(pp) + G.one * M, Pi = G.Slash(p) + G.one * M;
  const Mat Kin = helicity > 0 ? G.Slash(k) : G.Slash(kp) - G.one * mLep;   // spinor sum of "incoming" line
  const Mat Kout = helicity > 0 ? G.Slash(kp) + G.one * mLep : G.Slash(k);  // spinor sum of "outgoing" line
  cplx total = 0;
  for (int mu = 0; mu < 4; ++mu)
    for (int nu = 0; nu < 4; ++nu) {
      const cplx L = (Kout * G.g[mu] * PL * Kin * G.g[nu] * PL).Trace();
      const cplx H = (Pp * Gam[mu] * Pi * GamBar[nu]).Trace();
      total += L * H * Metric(mu) * Metric(nu);
    }
  return total.real();
}

}  // namespace nuflow::test
