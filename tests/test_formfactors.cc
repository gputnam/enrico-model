#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "core/Constants.h"
#include "kinematics/FormFactors.h"

using namespace nuflow;
using Catch::Approx;

TEST_CASE("Weak form factors at Q2 = 0", "[formfactors]") {
  for (const char* vec : {"BBBA05", "dipole"}) {
    YAML::Node n;
    n["vector"] = vec;
    n["axial"] = "dipole";
    n["MA"] = 1.03;
    n["gA"] = 1.2723;
    const auto ff = NucleonFormFactors::FromConfig(n).Evaluate(0.0);
    CHECK(ff.F1 == Approx(1.0).margin(1e-12));
    CHECK(ff.F2 == Approx(constants::kMuProton - constants::kMuNeutron - 1).epsilon(1e-12));
    CHECK(ff.FA == Approx(1.2723).epsilon(1e-12));
    // PCAC: F_P(0) = 2 M^2 g_A / m_pi^2
    const double M = constants::kMassNucleon, mpi = constants::kMassPionCharged;
    CHECK(ff.FP == Approx(2 * M * M * 1.2723 / (mpi * mpi)).epsilon(1e-12));
  }
}

TEST_CASE("BBBA05 Sachs form factors", "[formfactors]") {
  BBBA05FormFactors bbba;
  const auto s0 = bbba.Sachs(0);
  CHECK(s0.GEp == Approx(1.0));
  CHECK(s0.GMp == Approx(constants::kMuProton));
  CHECK(s0.GEn == Approx(0.0).margin(1e-15));
  CHECK(s0.GMn == Approx(constants::kMuNeutron));
  // Close to the dipole at moderate Q2 (within ~10%), and G_En small and positive.
  for (double Q2 : {0.1, 0.5, 1.0}) {
    const auto s = bbba.Sachs(Q2);
    const double GD = 1.0 / std::pow(1 + Q2 / 0.71, 2);
    CHECK(std::abs(s.GMp / (constants::kMuProton * GD) - 1) < 0.1);
    CHECK(s.GEn > 0);
    CHECK(s.GEn < 0.1);
  }
}

TEST_CASE("z-expansion axial form factor", "[formfactors]") {
  const double gA = 1.2723, mpi = 0.14;
  ZExpansionAxialFormFactor z(gA, {-2.30, 0.6, 3.8, -2.3}, 9 * mpi * mpi, -0.28);
  const auto& a = z.Coefficients();
  CHECK(z.FA(0) == Approx(gA).epsilon(1e-12));
  // Sum rules sum_k k^n a_k = 0, n = 0..3
  for (int n = 0; n < 4; ++n) {
    double s = 0;
    for (int k = 0; k <= 8; ++k) s += std::pow(double(k), n) * a[k];
    CHECK(s == Approx(0.0).margin(1e-9));
  }
  // The published central values (Meyer et al. Eq. 31, sign flipped to F_A(0) > 0) are rounded, so
  // individual high-order coefficients are not comparable (the sum rules amplify the rounding of
  // a1..a4); the form factor itself must agree.
  const double pub[9] = {0.759, -2.30, 0.6, 3.8, -2.3, -2.16, 0.896, 1.58, -0.823};
  const double s0 = std::sqrt(9 * mpi * mpi + 0.28);
  CHECK(a[0] == Approx(0.759).margin(0.01));
  for (double Q2 = 0.0; Q2 <= 1.0; Q2 += 0.1) {
    const double r = std::sqrt(9 * mpi * mpi + Q2);
    const double zz = (r - s0) / (r + s0);
    double fpub = 0;
    for (int k = 8; k >= 0; --k) fpub = fpub * zz + pub[k];
    INFO("Q2 = " << Q2);
    CHECK(z.FA(Q2) == Approx(fpub).epsilon(0.03));
  }
  // Falls with Q2 and stays within ~15% of a dipole with M_A = 1.0 GeV for Q2 < 1 GeV^2.
  double prev = z.FA(0);
  for (double Q2 = 0.05; Q2 <= 1.0; Q2 += 0.05) {
    const double f = z.FA(Q2);
    CHECK(f < prev);
    prev = f;
    CHECK(std::abs(f / (gA / std::pow(1 + Q2, 2)) - 1) < 0.15);
  }
}
