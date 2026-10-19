#include "kinematics/QECurrent.h"

#include "core/Constants.h"

namespace nuflow {

double ContractQE(const FourVector& k, const FourVector& kp, const FourVector& p, const FourVector& pp,
                  const WeakFormFactors& ff, double M, int helicity) {
  const FourVector q = pp - p;
  const FourVector P = pp + p;
  const double Q2 = -q.M2();
  const double M2 = M * M;

  // Vertex Gamma^mu = a gamma^mu + b P^mu + c gamma^mu gamma5 + d q^mu gamma5 (Gordon-decomposed).
  const double a = ff.F1 + ff.F2;
  const double b = -ff.F2 / (2 * M);
  const double c = ff.FA;
  const double d = ff.FP / M;

  const double kk = k.Dot(kp);
  // Symmetric lepton tensor contracted with  V^mu W^nu + W^mu V^nu - g^{mu nu} Y
  auto S = [&](const FourVector& V, const FourVector& W, double Y) {
    return 16 * (kp.Dot(V) * k.Dot(W) + kp.Dot(W) * k.Dot(V) - kk * V.Dot(W) + kk * Y);
  };
  // ... and with  Z^mu Z^nu
  auto R = [&](const FourVector& Z) { return 8 * (2 * kp.Dot(Z) * k.Dot(Z) - kk * Z.M2()); };

  const double sym = 4 * a * a * S(pp, p, 0.5 * Q2) + 4 * c * c * S(pp, p, 2 * M2 + 0.5 * Q2) +
                     (8 * a * b * M + 4 * b * b * (2 * M2 + 0.5 * Q2)) * R(P) + (2 * Q2 * d * d - 8 * M * c * d) * R(q);
  // Vector-axial interference (antisymmetric parts), sign fixed by helicity.
  const double anti = 128 * helicity * a * c * (k.Dot(p) * kp.Dot(q) - k.Dot(q) * kp.Dot(p));
  return sym + anti;
}

double LlewellynSmithDSigmaDQ2(double enu, double Q2, double m, double M, const WeakFormFactors& ff,
                               int helicity) {
  using namespace constants;
  const double M2 = M * M, m2 = m * m;
  const double tau = Q2 / (4 * M2);
  const double F1 = ff.F1, F2 = ff.F2, FA = ff.FA, FP = ff.FP;
  const double A = (m2 + Q2) / M2 *
                   ((1 + tau) * FA * FA - (1 - tau) * F1 * F1 + tau * (1 - tau) * F2 * F2 + 4 * tau * F1 * F2 -
                    m2 / (4 * M2) * ((F1 + F2) * (F1 + F2) + (FA + 2 * FP) * (FA + 2 * FP) - (Q2 / M2 + 4) * FP * FP));
  const double B = Q2 / M2 * FA * (F1 + F2);
  const double C = 0.25 * (FA * FA + F1 * F1 + tau * F2 * F2);
  const double su = 4 * M * enu - Q2 - m2;
  return kGF * kGF * kVud * kVud * M2 / (8 * kPi * enu * enu) *
         (A + helicity * B * su / M2 + C * su * su / (M2 * M2));
}

}  // namespace nuflow
