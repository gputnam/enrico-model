#include "fsi/HadronNucleonXSec.h"

#include <algorithm>
#include <cmath>

#include "core/Constants.h"

namespace nuflow::hadxs {

using constants::kMassPionAverage;
using constants::kMassNucleon;
using constants::kPi;

// Cugnon, L'Hote, Vandermeulen (1996). p in GeV/c, sigma in mb. Below p = 0.1 GeV/c the
// parametrisation is frozen (its stated range of validity); such collisions are Pauli blocked anyway.
double NNTotal(bool like, double p) {
  p = std::max(p, 0.1);
  if (like) {
    if (p < 0.44) return 34.0 * std::pow(p / 0.4, -2.104);
    if (p < 0.8) return 23.5 + 1000.0 * std::pow(p - 0.7, 4);
    if (p < 1.5) return 23.5 + 24.6 / (1.0 + std::exp(-(p - 1.2) / 0.1));
    if (p < 5.0) return 41.0 + 60.0 * (p - 0.9) * std::exp(-1.2 * p);
    const double l = std::log(p);
    return 48.0 + 0.522 * l * l - 4.51 * l;
  }
  if (p < 0.44) return 6.3555 * std::pow(p, -3.2481) * std::exp(-0.377 * std::log(p) * std::log(p));
  if (p < 1.0) return 33.0 + 196.0 * std::pow(std::abs(p - 0.95), 2.5);
  if (p < 2.0) return 24.2 + 8.9 * p;
  return 42.0;
}

double NNElastic(bool like, double p) {
  p = std::max(p, 0.1);
  double el;
  if (like) {
    if (p < 0.44) el = 34.0 * std::pow(p / 0.4, -2.104);
    else if (p < 0.8) el = 23.5 + 1000.0 * std::pow(p - 0.7, 4);
    else if (p < 2.0) el = 1250.0 / (p + 50.0) - 4.0 * (p - 1.3) * (p - 1.3);
    else el = 77.0 / (p + 1.5);
  } else {
    if (p < 0.44) el = 6.3555 * std::pow(p, -3.2481) * std::exp(-0.377 * std::log(p) * std::log(p));
    else if (p < 0.8) el = 33.0 + 196.0 * std::pow(std::abs(p - 0.95), 2.5);
    else if (p < 2.0) el = 31.0 / std::sqrt(p);
    else el = 77.0 / (p + 1.5);
  }
  return std::min(el, NNTotal(like, p));
}

double NNSlope(bool like, double p) {
  auto ppSlope = [](double q) {
    if (q < 2.0) {
      const double q8 = std::pow(q, 8);
      return 5.5 * q8 / (7.7 + q8);
    }
    return 5.334 + 0.67 * (q - 2.0);
  };
  if (like) return ppSlope(p);
  if (p < 0.225) return 0.0;
  if (p < 0.6) return 16.53 * (p - 0.225);
  if (p < 1.6) return -1.63 * p + 7.16;
  return ppSlope(p);
}

double CMMomentum(double W, double m1, double m2) {
  const double s = W * W;
  const double lam = (s - (m1 + m2) * (m1 + m2)) * (s - (m1 - m2) * (m1 - m2));
  return lam > 0 ? std::sqrt(lam) / (2 * W) : 0.0;
}

double DeltaWidth(double W) {
  // HNV / CCRES: Gamma = (1/6pi) (f*/m_pi)^2 (M/W) q_cm^3 (theory/pion_production.md).
  const double q = CMMomentum(W, kMassNucleon, kMassPionAverage);
  if (q <= 0) return 0;
  const double f = kFStar / kMassPionAverage;
  return f * f / (6 * kPi) * (kMassNucleon / W) * q * q * q;
}

double DeltaRegionWeight(double W) {
  if (W <= 1.35) return 1.0;
  if (W >= 1.45) return 0.0;
  return (1.45 - W) / 0.10;
}

double DeltaIsospinFactor(int qpi, int qN) {
  // |<3/2, qpi+qN-1/2 | 1 qpi ; 1/2 qN-1/2>|^2
  const int qD = qpi + qN;
  if (qD == 2 || qD == -1) return 1.0;      // pi+ p, pi- n
  if (qpi == 0) return 2.0 / 3.0;          // pi0 p -> Delta+, pi0 n -> Delta0
  return 1.0 / 3.0;                        // pi+ n -> Delta+, pi- p -> Delta0
}

OsetSelfEnergy DeltaSelfEnergy(double tpi, double rr) {
  // Oset & Salcedo, NPA 468 (1987) 631: C(x) = a x^2 + b x + c [MeV], x = T_pi / m_pi.
  const double mpi = constants::kMassPionCharged;
  const double x = std::clamp(tpi, 0.085, 0.315) / mpi;
  auto poly = [x](double a, double b, double c) { return a * x * x + b * x + c; };
  const double cQ = poly(-5.19, 15.35, 2.06), cA2 = poly(1.06, -6.64, 22.66), cA3 = poly(-13.46, 46.17, -20.34);
  const double alpha = poly(0.382, -1.322, 1.466), beta = poly(-0.038, 0.204, 0.613);
  OsetSelfEnergy s;
  if (rr <= 0) return s;
  s.cQ = 1e-3 * std::max(0.0, cQ) * std::pow(rr, alpha);
  s.cA2 = 1e-3 * std::max(0.0, cA2) * std::pow(rr, beta);
  s.cA3 = 1e-3 * std::max(0.0, cA3) * std::pow(rr, 2 * beta);
  return s;
}

}  // namespace nuflow::hadxs
