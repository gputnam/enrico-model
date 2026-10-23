#include "kinematics/FormFactors.h"

#include <cmath>
#include <stdexcept>

#include "core/Config.h"
#include "core/Constants.h"

namespace nuflow {

using namespace constants;

namespace {

double Dipole(double Q2, double M2) {
  const double d = 1.0 + Q2 / M2;
  return 1.0 / (d * d);
}

// BBBA05 rational form: (sum_{k=0}^{2} a_k tau^k) / (1 + sum_{k=1}^{4} b_k tau^k)
double Kelly(double tau, const double a[3], const double b[4]) {
  const double num = a[0] + tau * (a[1] + tau * a[2]);
  const double den = 1.0 + tau * (b[0] + tau * (b[1] + tau * (b[2] + tau * b[3])));
  return num / den;
}

}  // namespace

SachsFormFactors DipoleGalsterFormFactors::Sachs(double Q2) const {
  const double GD = Dipole(Q2, MV2_);
  const double tauN = Q2 / (4 * kMassNeutron * kMassNeutron);
  SachsFormFactors s;
  s.GEp = GD;
  s.GMp = kMuProton * GD;
  s.GMn = kMuNeutron * GD;
  s.GEn = -kMuNeutron * tauN / (1.0 + 5.6 * tauN) * GD;
  return s;
}

SachsFormFactors BBBA05FormFactors::Sachs(double Q2) const {
  // Table 1 of hep-ex/0602017; tau = Q2 / 4M^2 with a single nucleon mass M.
  static const double aEp[3] = {1.0, -0.0578, 0.0}, bEp[4] = {11.1, 13.6, 33.0, 0.0};
  static const double aMp[3] = {1.0, 0.150, 0.0}, bMp[4] = {11.1, 19.6, 7.54, 0.0};
  static const double aEn[3] = {0.0, 1.25, 1.30}, bEn[4] = {-9.86, 305.0, -758.0, 802.0};
  static const double aMn[3] = {1.0, 1.81, 0.0}, bMn[4] = {14.1, 20.7, 68.7, 0.0};
  const double tau = Q2 / (4 * kMassNucleon * kMassNucleon);
  SachsFormFactors s;
  s.GEp = Kelly(tau, aEp, bEp);
  s.GMp = kMuProton * Kelly(tau, aMp, bMp);
  s.GEn = Kelly(tau, aEn, bEn);
  s.GMn = kMuNeutron * Kelly(tau, aMn, bMn);
  return s;
}

ZExpansionAxialFormFactor::ZExpansionAxialFormFactor(double gA, const std::vector<double>& aFree, double tcut,
                                                     double t0)
    : tcut_(tcut), t0_(t0) {
  if (aFree.empty()) throw std::runtime_error("z-expansion: need at least one free coefficient a_1");
  if (tcut <= 0 || t0 >= tcut) throw std::runtime_error("z-expansion: need tcut > 0 and t0 < tcut");
  const int nFree = int(aFree.size());
  const int kmax = nFree + 4;
  a_.assign(kmax + 1, 0.0);
  // Unknowns x = (a0, a_{kmax-3}, ..., a_kmax). Row 0: F_A(0) = gA. Rows 1..4: sum_k k^n a_k = 0, n = 0..3.
  const double z0 = Z(0.0);
  const int unknown[5] = {0, kmax - 3, kmax - 2, kmax - 1, kmax};
  double m[5][6] = {};
  for (int j = 0; j < 5; ++j) m[0][j] = std::pow(z0, unknown[j]);
  double rhs0 = gA;
  for (int k = 1; k <= nFree; ++k) rhs0 -= aFree[k - 1] * std::pow(z0, k);
  m[0][5] = rhs0;
  for (int n = 0; n < 4; ++n) {
    for (int j = 0; j < 5; ++j) m[n + 1][j] = std::pow(double(unknown[j]), n);  // 0^0 = 1
    double r = 0;
    for (int k = 1; k <= nFree; ++k) r -= std::pow(double(k), n) * aFree[k - 1];
    m[n + 1][5] = r;
  }
  // Gaussian elimination with partial pivoting.
  for (int c = 0; c < 5; ++c) {
    int piv = c;
    for (int r = c + 1; r < 5; ++r)
      if (std::abs(m[r][c]) > std::abs(m[piv][c])) piv = r;
    for (int j = 0; j < 6; ++j) std::swap(m[c][j], m[piv][j]);
    if (std::abs(m[c][c]) < 1e-14) throw std::runtime_error("z-expansion: singular sum-rule system");
    for (int r = 0; r < 5; ++r) {
      if (r == c) continue;
      const double f = m[r][c] / m[c][c];
      for (int j = c; j < 6; ++j) m[r][j] -= f * m[c][j];
    }
  }
  for (int k = 1; k <= nFree; ++k) a_[k] = aFree[k - 1];
  for (int j = 0; j < 5; ++j) a_[unknown[j]] = m[j][5] / m[j][j];
}

double ZExpansionAxialFormFactor::Z(double Q2) const {
  const double s1 = std::sqrt(tcut_ + Q2), s0 = std::sqrt(tcut_ - t0_);
  return (s1 - s0) / (s1 + s0);
}

double ZExpansionAxialFormFactor::FA(double Q2) const {
  const double z = Z(Q2);
  double r = 0;
  for (int k = int(a_.size()) - 1; k >= 0; --k) r = r * z + a_[k];
  return r;
}

WeakFormFactors NucleonFormFactors::Evaluate(double Q2) const {
  const SachsFormFactors s = vector_->Sachs(Q2);
  // CVC: isovector combinations of the electromagnetic form factors.
  const double GEV = s.GEp - s.GEn;
  const double GMV = s.GMp - s.GMn;
  const double tau = Q2 / (4 * kMassNucleon * kMassNucleon);
  WeakFormFactors f;
  f.F1 = (GEV + tau * GMV) / (1 + tau);
  f.F2 = (GMV - GEV) / (1 + tau);
  f.FA = axial_->FA(Q2);
  // PCAC / pion-pole dominance.
  f.FP = 2 * kMassNucleon * kMassNucleon * f.FA / (kMassPionCharged * kMassPionCharged + Q2);
  return f;
}

NucleonFormFactors NucleonFormFactors::FromConfig(const YAML::Node& node) {
  std::shared_ptr<const VectorFormFactorModel> vec;
  const std::string v = Get<std::string>(node, "vector", "BBBA05");
  if (v == "BBBA05") vec = std::make_shared<BBBA05FormFactors>();
  else if (v == "dipole") vec = std::make_shared<DipoleGalsterFormFactors>(Get<double>(node, "MV2", 0.71));
  else throw std::runtime_error("form factors: unknown vector model '" + v + "'");

  const double gA = Get<double>(node, "gA", 1.2723);
  std::shared_ptr<const AxialFormFactorModel> ax;
  const std::string a = Get<std::string>(node, "axial", "dipole");
  if (a == "dipole") {
    ax = std::make_shared<DipoleAxialFormFactor>(gA, Get<double>(node, "MA", 1.03));
  } else if (a == "zexp") {
    // Defaults: central values of the deuterium fit of Meyer et al. (sign flipped to F_A(0) > 0).
    // zexp_a lists the free coefficients a_1..a_{kmax-4}; kmax = length + 4 (4 sum rules).
    std::vector<double> c = {-2.30, 0.6, 3.8, -2.3};
    if (node["zexp_a"]) {
      c = node["zexp_a"].as<std::vector<double>>();
      if (c.empty()) throw std::runtime_error("form factors: zexp_a needs at least a1");
    }
    const double mpi = 0.14;
    ax = std::make_shared<ZExpansionAxialFormFactor>(gA, c, Get<double>(node, "tcut", 9 * mpi * mpi),
                                                     Get<double>(node, "t0", -0.28));
  } else {
    throw std::runtime_error("form factors: unknown axial model '" + a + "'");
  }
  return NucleonFormFactors(vec, ax);
}

}  // namespace nuflow
