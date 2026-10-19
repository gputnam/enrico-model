#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "core/Nucleus.h"
#include "core/PDG.h"
#include "initial_state/FreeNucleon.h"
#include "initial_state/RelativisticFermiGas.h"

using namespace nuflow;
using Catch::Approx;

TEST_CASE("Nuclear database", "[nuclear]") {
  const auto& db = NuclearDatabase::Instance();
  const Nucleus& C = db.Get(1000060120);
  CHECK(C.Z == 6);
  CHECK(C.A == 12);
  CHECK(C.kFermiNeutron == Approx(0.221));
  CHECK(C.removalEnergy == Approx(0.025));
  // Hydrogen is accepted both as 2212 and as 1000010010.
  CHECK(db.Get(1000010010).pdg == pdg::kProton);
  CHECK(db.Get(pdg::kProton).IsFreeNucleon());
}

TEST_CASE("Relativistic Fermi gas momentum distribution", "[nuclear]") {
  const Nucleus& C = NuclearDatabase::Instance().Get(1000060120);
  RelativisticFermiGas rfg;
  Random rng(3);
  HitNucleonRequest req;
  req.species[0] = pdg::kNeutron;
  const int N = 200000;
  double sumK2 = 0, maxK = 0;
  for (int i = 0; i < N; ++i) {
    InitialState s;
    rfg.Sample(C, req, s, rng);
    const double k = s.hit[0].p.Mag();
    sumK2 += k * k;
    maxK = std::max(maxK, k);
    // Smith-Moniz energy and exact four-momentum balance with the remnant
    CHECK(s.hit[0].energy == Approx(std::sqrt(k * k + pdg::NucleonMass(pdg::kNeutron) * pdg::NucleonMass(pdg::kNeutron)) - C.removalEnergy));
    const FourVector tot = s.hit[0].Momentum() + s.remnant;
    REQUIRE(std::abs(tot.e - C.mass) < 1e-9);
    REQUIRE(tot.P() < 1e-12);
    REQUIRE(s.remnantPdg == pdg::IonCode(6, 11));
  }
  const double kF = C.kFermiNeutron;
  CHECK(maxK <= kF);
  CHECK(sumK2 / N == Approx(0.6 * kF * kF).epsilon(0.01));  // <k^2> = 3/5 kF^2
  CHECK(rfg.PauliBlocked(C, pdg::kProton, {0, 0, 0.9 * kF}));
  CHECK_FALSE(rfg.PauliBlocked(C, pdg::kProton, {0, 0, 1.1 * kF}));
}

TEST_CASE("Fermi gas with high-momentum tail", "[nuclear]") {
  const Nucleus& C = NuclearDatabase::Instance().Get(1000060120);
  RelativisticFermiGas::Options o;
  o.tailFraction = 0.2;
  RelativisticFermiGas rfg(o);
  Random rng(5);
  int above = 0;
  const int N = 100000;
  for (int i = 0; i < N; ++i) {
    const double k = rfg.SampleMomentum(C.kFermiNeutron, rng);
    REQUIRE(k <= o.tailCutoff);
    above += k > C.kFermiNeutron;
  }
  CHECK(double(above) / N == Approx(0.2).margin(0.005));

  // Bodek-Ritchie: tail fraction 6 (kF a / pi)^2 = 11.9% for carbon, 1/k^4 up to 4 GeV.
  RelativisticFermiGas::Options br;
  br.bodekRitchie = true;
  RelativisticFermiGas brModel(br);
  CHECK(brModel.TailFraction(0.221) == Approx(6 * std::pow(0.221 * 2 / M_PI, 2)));
  CHECK(brModel.TailFraction(0.221) == Approx(0.119).margin(0.001));
  above = 0;
  int above1 = 0;
  for (int i = 0; i < N; ++i) {
    const double k = brModel.SampleMomentum(0.221, rng);
    REQUIRE(k <= 4.0);
    above += k > 0.221;
    above1 += k > 1.0;
  }
  CHECK(double(above) / N == Approx(0.119).margin(0.004));
  // P(k > 1 | tail) = (1/1 - 1/4) / (1/kF - 1/4)
  CHECK(double(above1) / above == Approx((1.0 - 0.25) / (1 / 0.221 - 0.25)).margin(0.01));
}
