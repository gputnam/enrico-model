#pragma once
// Explicit 4x4 complex Dirac matrices (Dirac representation), for matrix elements that are evaluated
// numerically as traces rather than in closed form (e.g. pion production, theory/pion_production.md).
//
// Conventions (theory/README.md): metric (+,-,-,-), gamma5 = i gamma^0 gamma^1 gamma^2 gamma^3,
// sigma^{mu nu} = (i/2)[gamma^mu, gamma^nu], Bar(X) = gamma^0 X^dagger gamma^0.

#include <array>
#include <complex>

#include "core/FourVector.h"

namespace nuflow::dirac {

using cplx = std::complex<double>;

struct Mat {
  std::array<cplx, 16> m{};

  cplx& operator()(int i, int j) { return m[4 * i + j]; }
  const cplx& operator()(int i, int j) const { return m[4 * i + j]; }

  Mat operator*(const Mat& o) const {
    Mat r;
    for (int i = 0; i < 4; ++i)
      for (int k = 0; k < 4; ++k) {
        const cplx a = m[4 * i + k];
        if (a == cplx(0)) continue;
        for (int j = 0; j < 4; ++j) r.m[4 * i + j] += a * o.m[4 * k + j];
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
  Mat& operator+=(const Mat& o) {
    for (int i = 0; i < 16; ++i) m[i] += o.m[i];
    return *this;
  }
  // this += s * o
  void AddScaled(const Mat& o, cplx s) {
    for (int i = 0; i < 16; ++i) m[i] += s * o.m[i];
  }
  Mat Dagger() const {
    Mat r;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) r.m[4 * i + j] = std::conj(m[4 * j + i]);
    return r;
  }
  cplx Trace() const { return m[0] + m[5] + m[10] + m[15]; }
  static Mat Identity(cplx s = 1.0) {
    Mat r;
    for (int i = 0; i < 4; ++i) r.m[5 * i] = s;
    return r;
  }
};

// Tr[A B] without forming the product.
inline cplx TraceProduct(const Mat& a, const Mat& b) {
  cplx s = 0;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) s += a.m[4 * i + j] * b.m[4 * j + i];
  return s;
}

inline double Metric(int mu) { return mu == 0 ? 1.0 : -1.0; }
inline double Comp(const FourVector& v, int mu) { return mu == 0 ? v.e : mu == 1 ? v.px : mu == 2 ? v.py : v.pz; }

struct Gammas {
  Mat g[4], g5, one, g5g[4];  // g5g[mu] = gamma^mu gamma5
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
    for (int mu = 0; mu < 4; ++mu) g5g[mu] = g[mu] * g5;
  }
  Mat Slash(const FourVector& v) const {
    Mat r;
    r.AddScaled(g[0], v.e);
    r.AddScaled(g[1], -v.px);
    r.AddScaled(g[2], -v.py);
    r.AddScaled(g[3], -v.pz);
    return r;
  }
  Mat Bar(const Mat& x) const { return g[0] * x.Dagger() * g[0]; }

  static const Gammas& Get() {
    static const Gammas G;
    return G;
  }
};

}  // namespace nuflow::dirac
