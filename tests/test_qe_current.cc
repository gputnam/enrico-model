// Validation of the CCQE matrix element:
//  1. closed-form contraction == explicit Dirac traces (Standard-Model V-A current), nu and nubar;
//  2. free nucleon at rest: closed form == Llewellyn Smith dsigma/dQ2;
//  3. generator phase space: Monte Carlo total cross section on a free nucleon == integral of
//     the Llewellyn Smith dsigma/dQ2 (tests the CM sampling, Jacobian, flux factor and units).
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "DiracAlgebra.h"
#include "core/Constants.h"
#include "core/Nucleus.h"
#include "core/PDG.h"
#include "core/Random.h"
#include "initial_state/FreeNucleon.h"
#include "kinematics/QECurrent.h"
#include "kinematics/QuasiElasticKinematics.h"

using namespace nuflow;
using Catch::Approx;

namespace {
FourVector OnShell(double M, const ThreeVector& p) { return {std::sqrt(p.Mag2() + M * M), p}; }
ThreeVector RandomVec(Random& r, double scale) {
  return {scale * (2 * r.Uniform() - 1), scale * (2 * r.Uniform() - 1), scale * (2 * r.Uniform() - 1)};
}
}  // namespace

TEST_CASE("Closed-form QE contraction matches explicit Dirac traces", "[qe]") {
  Random rng(42);
  const double M = constants::kMassNucleon, m = constants::kMassMuon;
  for (int trial = 0; trial < 50; ++trial) {
    const WeakFormFactors ff{rng.Uniform(-1, 2), rng.Uniform(-1, 4), rng.Uniform(-2, 2), rng.Uniform(-5, 20)};
    const FourVector k = OnShell(0, RandomVec(rng, 3));
    const FourVector kp = OnShell(m, RandomVec(rng, 2));
    const FourVector p = OnShell(M, RandomVec(rng, 0.3));
    const FourVector pp = OnShell(M, RandomVec(rng, 1.5));
    for (int h : {+1, -1}) {
      const double closed = ContractQE(k, kp, p, pp, ff, M, h);
      const double brute = test::BruteForceContraction(k, kp, m, p, pp, M, ff, h, /*axialSign=*/-1);
      CHECK(closed == Approx(brute).epsilon(1e-10));
    }
  }
}

TEST_CASE("Free nucleon at rest reproduces Llewellyn Smith", "[qe]") {
  const double M = constants::kMassNucleon, m = constants::kMassMuon;
  for (int h : {+1, -1})
    for (double E : {0.3, 1.0, 5.0})
      for (double Q2 : {0.05, 0.3, 0.8}) {
        const double El = E - Q2 / (2 * M);
        if (El <= m) continue;
        const double pl = std::sqrt(El * El - m * m);
        const double cth = (2 * E * El - m * m - Q2) / (2 * E * pl);
        if (std::abs(cth) > 1) continue;
        const FourVector k(E, 0, 0, E), p(M, 0, 0, 0);
        const FourVector kp(El, pl * std::sqrt(1 - cth * cth), 0, pl * cth);
        const WeakFormFactors ff{0.8, 2.5, 1.0, 3.0};
        const double LH = ContractQE(k, kp, p, k + p - kp, ff, M, h);
        // dsigma/dQ2 = <|M|^2> / (64 pi M^2 E^2),  <|M|^2> = G_F^2 Vud^2 LH / 4
        const double ds = constants::kGF * constants::kGF * constants::kVud * constants::kVud * LH /
                          (256 * constants::kPi * M * M * E * E);
        CHECK(ds == Approx(LlewellynSmithDSigmaDQ2(E, Q2, m, M, ff, h)).epsilon(1e-10));
      }
}

namespace {
// Integral of the Llewellyn Smith formula over the kinematically allowed Q2 range for a free
// nucleon at rest (equal masses M), Simpson's rule.
double IntegratedLS(double E, const NucleonFormFactors& ffm, int h) {
  const double M = constants::kMassNucleon, m = constants::kMassMuon;
  const double s = M * M + 2 * M * E, rs = std::sqrt(s);
  const double kcm = (s - M * M) / (2 * rs);
  const double El = (s + m * m - M * M) / (2 * rs), pl = std::sqrt(El * El - m * m);
  const double q2min = -m * m + 2 * kcm * (El - pl), q2max = -m * m + 2 * kcm * (El + pl);
  const int n = 4000;
  double sum = 0;
  for (int i = 0; i <= n; ++i) {
    const double Q2 = q2min + (q2max - q2min) * i / n;
    const double w = (i == 0 || i == n) ? 1 : (i % 2 ? 4 : 2);
    sum += w * LlewellynSmithDSigmaDQ2(E, Q2, m, M, ffm.Evaluate(Q2), h);
  }
  return sum * (q2max - q2min) / (3 * n) * constants::kGeV2To1e38Cm2;
}
}  // namespace

TEST_CASE("Monte Carlo free-nucleon cross section equals integrated Llewellyn Smith", "[qe]") {
  YAML::Node cfg;
  cfg["vector"] = "BBBA05";
  cfg["axial"] = "dipole";
  cfg["MA"] = 1.03;
  const NucleonFormFactors ffm = NucleonFormFactors::FromConfig(cfg);
  QuasiElasticKinematics kin(ffm);
  FreeNucleon model;
  const auto& db = NuclearDatabase::Instance();
  for (int nu : {14, -14}) {
    const Nucleus& target = db.Get(nu > 0 ? pdg::kNeutron : pdg::kProton);
    HitNucleonRequest req;
    req.species[0] = nu > 0 ? pdg::kNeutron : pdg::kProton;
    for (double E : {0.4, 1.0, 3.0}) {
      Random rng(7, static_cast<uint64_t>(E * 1000) + (nu > 0));
      Event ev;
      InitialState init;
      model.Sample(target, req, init, rng);
      const long N = 400000;
      double sum = 0, sum2 = 0;
      for (long i = 0; i < N; ++i) {
        const double w = kin.Generate(nu, E, target, init, model, ev, rng);
        sum += w;
        sum2 += w * w;
      }
      const double mean = sum / N, err = std::sqrt((sum2 / N - mean * mean) / N);
      const double ref = IntegratedLS(E, ffm, nu > 0 ? 1 : -1);
      INFO("nu=" << nu << " E=" << E << " MC=" << mean << " +- " << err << " LS=" << ref);
      // The generator uses the physical n and p masses while the LS closed form uses their average;
      // this shifts nu-n by up to ~0.6% (with equal masses the two agree to < 0.05%, see
      // theory/quasielastic.md, validation section).
      CHECK(std::abs(mean - ref) < 3 * err + 7e-3 * ref);
    }
  }
}
