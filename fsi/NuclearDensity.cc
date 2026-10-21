#include "fsi/NuclearDensity.h"

#include <cmath>
#include <stdexcept>

#include "core/Constants.h"
#include "core/PDG.h"

namespace nuflow {

namespace {
constexpr double kHbarC = 0.1973269804;      // GeV fm
constexpr double kRp2 = 0.8409 * 0.8409;     // proton charge radius^2 [fm^2] (CODATA 2018)
constexpr double kRn2 = -0.1161;             // neutron mean-square charge radius [fm^2] (PDG)
// Darwin-Foldy term 3/(4 M^2) [fm^2]
constexpr double kDarwinFoldy = 0.75 * kHbarC * kHbarC / (constants::kMassProton * constants::kMassProton);
}  // namespace

NuclearDensity::NuclearDensity(const Nucleus& nucleus) : Z_(nucleus.Z), A_(nucleus.A) {
  if (nucleus.A < 2) throw std::runtime_error("NuclearDensity: not a nucleus");
  double rch = nucleus.rmsChargeRadius;
  if (rch <= 0) rch = 0.94 * std::cbrt(static_cast<double>(nucleus.A)) + 0.4;  // rough systematics
  const int N = nucleus.A - nucleus.Z;
  const double r2pt = rch * rch - kRp2 - (nucleus.Z > 0 ? double(N) / nucleus.Z : 0.0) * kRn2 - kDarwinFoldy;
  rmsPoint_ = std::sqrt(r2pt);
  ho_ = nucleus.A <= 16;
  if (ho_) {
    c_ = std::max(0.0, (nucleus.Z - 2) / 3.0);
    // <x^2> = (6 + 15 c) / (4 + 6 c)
    b_ = std::sqrt(r2pt * (4 + 6 * c_) / (6 + 15 * c_));
  } else {
    // 2pF: <r^2> ~ (3/5) R^2 + (7/5) pi^2 a^2 (leading terms)
    const double R2 = (r2pt - 1.4 * constants::kPi * constants::kPi * a_ * a_) / 0.6;
    b_ = std::sqrt(std::max(R2, 1.0));
  }
  // Normalise to A numerically and find rMax and the maximum of r^2 rho.
  double sum = 0;
  const int n = 20000;
  const double rTop = ho_ ? 8 * b_ : b_ + 20 * a_;
  const double dr = rTop / n;
  for (int i = 0; i < n; ++i) {
    const double r = (i + 0.5) * dr;
    sum += 4 * constants::kPi * r * r * Shape(r) * dr;
  }
  norm_ = nucleus.A / sum;
  const double rho0 = Rho(0);
  rMax_ = rTop;
  for (int i = 0; i < n; ++i) {
    const double r = (i + 0.5) * dr;
    maxR2Rho_ = std::max(maxR2Rho_, r * r * Rho(r));
    if (Rho(r) < 1e-4 * rho0) {
      rMax_ = r;
      break;
    }
  }
  maxR2Rho_ *= 1.01;
}

double NuclearDensity::Shape(double r) const {
  if (ho_) {
    const double x2 = r * r / (b_ * b_);
    return (1 + c_ * x2) * std::exp(-x2);
  }
  return 1.0 / (1.0 + std::exp((r - b_) / a_));
}

double NuclearDensity::Rho(double r) const { return norm_ * Shape(r); }

double NuclearDensity::RhoSpecies(int pdgCode, double r) const {
  const double f = pdgCode == pdg::kProton ? double(Z_) / A_ : double(A_ - Z_) / A_;
  return f * Rho(r);
}

double NuclearDensity::KFermi(int pdgCode, double r) const {
  return kHbarC * std::cbrt(3 * constants::kPi * constants::kPi * RhoSpecies(pdgCode, r));
}

ThreeVector NuclearDensity::SamplePosition(Random& rng) const {
  double r;
  do {
    r = rMax_ * rng.Uniform();
  } while (rng.Uniform() * maxR2Rho_ > r * r * Rho(r));
  const double cosT = 2 * rng.Uniform() - 1, sinT = std::sqrt(std::max(0.0, 1 - cosT * cosT));
  const double phi = 2 * constants::kPi * rng.Uniform();
  return {r * sinT * std::cos(phi), r * sinT * std::sin(phi), r * cosT};
}

}  // namespace nuflow
