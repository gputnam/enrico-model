// Validation of the CCCOH channel (Berger-Sehgal PCAC coherent pion production; theory/coherent_pion.md):
//  1. the pi-12C elastic input: BS09 table values, interpolation, A1/b1 continuity into the Rein-Sehgal
//     hadronic model, PDG pi N total cross section;
//  2. the lepton-mass factor reduces to G_A^2 for a massless lepton and is < G_A^2 at Q2 ~ Q2min;
//  3. events: four-momentum conservation with the recoiling ground-state nucleus, exactly one charged
//     lepton and one pion of the right charge as status-1 particles, no nucleons, no cascade, process 100,
//     and the channel does not apply to hydrogen.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "core/Constants.h"
#include "core/PDG.h"
#include "kinematics/CoherentPionKinematics.h"
#include "orchestration/Generator.h"

using namespace nuflow;
using Catch::Approx;

namespace {

struct MemWriter : OutputWriter {
  std::vector<Event> events;
  void Begin(const RunInfo&) override {}
  void Write(const Event& e) override { events.push_back(e); }
  void End(const RunSummary&) override {}
};

YAML::Node CohConfig(const std::string& target, double e, int nu, const std::string& fsi) {
  YAML::Node c = YAML::Load(R"(
run: {seed: 11, batch_size: 500, threads: 4}
channels:
  - type: CCCOH
xsec_table: {nodes: 2, samples: 100000, cache: none}
)");
  c["fsi"]["type"] = fsi;
  c["beam"]["pdg"] = nu;
  c["beam"]["flux"]["type"] = "mono";
  c["beam"]["flux"]["energy"] = e;
  c["target"] = target;
  return c;
}

}  // namespace

TEST_CASE("pi-12C elastic input of the coherent model", "[coherent]") {
  PionNucleusElastic el;
  double A1, b1;
  el.Coefficients(0.162, 12, A1, b1);
  CHECK(A1 == Approx(22400));
  CHECK(b1 == Approx(89.2));
  el.Coefficients(0.87, 12, A1, b1);
  CHECK(A1 == Approx(5140));
  // continuity at both ends of the interpolation region
  double aR, bR;
  el.ReinSehgal(1.7, 12, aR, bR);
  el.Coefficients(1.7 - 1e-9, 12, A1, b1);
  CHECK(A1 == Approx(aR).epsilon(1e-6));
  // BS09: the RS hadronic model and the data meet near p_pi ~ 1 GeV at sigma_el ~ 80 mb
  el.ReinSehgal(0.87, 12, aR, bR);
  CHECK(aR / bR > 60);
  CHECK(aR / bR < 110);
  CHECK(el.Sigma(0.87, 12) == Approx(5140 / 62.2));
  // PDG isoscalar pi N total: ~25 mb at p = 10 GeV, ~24 mb at 100 GeV
  CHECK(PionNucleusElastic::PiNTotal(10.0) == Approx(25.0).margin(1.5));
  CHECK(PionNucleusElastic::PiNTotal(100.0) == Approx(24.0).margin(1.5));
  CHECK(PionNucleusElastic::PiNElastic(10.0) == Approx(4.6).margin(1.0));
  // A scaling: sigma_el ~ A^(2/3) in the data region
  CHECK(el.Sigma(0.3, 40) / el.Sigma(0.3, 12) == Approx(std::pow(40.0 / 12, 2.0 / 3)).epsilon(1e-9));
}

TEST_CASE("coherent lepton-mass (pion-pole) factor", "[coherent]") {
  CoherentPionKinematics kin(PionNucleusElastic{}, CoherentPionKinematics::Options{});
  CoherentPionKinematics::Options o;
  o.leptonMassCorrection = false;
  CoherentPionKinematics kin0(PionNucleusElastic{}, o);
  // nu_e (m_e ~ 0): identical with and without the correction
  CHECK(kin.DSigma(12, 3.0, 0.05, 1.0, 0.01, 0.8, 12) == Approx(kin0.DSigma(12, 3.0, 0.05, 1.0, 0.01, 0.8, 12)).epsilon(1e-4));
  // nu_mu at Q2 = Q2min: suppressed (Adler screening); at Q2 >> m_pi^2 -> ~ G_A^2
  const double y = 0.5, E = 1.0, Q2min = constants::kMassMuon * constants::kMassMuon * y / (1 - y);
  CHECK(kin.DSigma(14, E, Q2min, y * E, 0.01, 0.3, 12) < 0.8 * kin0.DSigma(14, E, Q2min, y * E, 0.01, 0.3, 12));
  CHECK(kin.DSigma(14, 6.0, 0.5, 3.0, 0.01, 2.9, 12) == Approx(kin0.DSigma(14, 6.0, 0.5, 3.0, 0.01, 2.9, 12)).epsilon(0.05));
}

TEST_CASE("CCCOH events: kinematics, final state, conservation", "[coherent]") {
  for (int nu : {14, -14})
    for (const char* fsi : {"cascade", "none"}) {
      Generator gen(CohConfig("1000060120:1,1000010010:1", 6.0, nu, fsi));
      MemWriter w;
      gen.Run(2000, {&w});
      double meanT = 0;
      for (const auto& ev : w.events) {
        REQUIRE(ev.channel == 100);
        REQUIRE(ev.targetPdg == 1000060120);
        FourVector in, out, rem;
        int nLep = 0, nPi = 0, nOther = 0, nRem = 0;
        for (const auto& p : ev.particles) {
          if (p.status == kStatusBeam || p.status == kStatusTarget) in = in + p.p;
          if (p.status == kStatusPreFSI) ++nOther;
          if (p.status == kStatusRemnantPostFSI || (p.status == kStatusRemnant && std::string(fsi) == "none")) {
            rem = p.p;
            ++nRem;
            REQUIRE(p.pdg == 1000060120);
          }
          if (p.status != kStatusFinal) continue;
          out = out + p.p;
          if (p.pdg == pdg::ChargedLeptonFor(nu)) ++nLep;
          else if (p.pdg == (nu > 0 ? 211 : -211)) ++nPi;
          else ++nOther;
        }
        REQUIRE(nLep == 1);
        REQUIRE(nPi == 1);
        REQUIRE(nOther == 0);
        REQUIRE(nRem == 1);
        const FourVector d = in - out - rem;
        REQUIRE(std::abs(d.e) < 1e-9);
        REQUIRE(d.P() < 1e-9);
        REQUIRE(rem.M() == Approx(gen.Info().target.components[0].nucleus.mass).epsilon(1e-10));
        const double t = -(rem - FourVector(rem.M(), 0, 0, 0)).M2();
        meanT += t;
        REQUIRE(t >= 0);
        REQUIRE(ev.Q2 >= 0);
      }
      meanT /= w.events.size();
      CHECK(meanT < 0.05);  // coherent: |t| ~ 1/b ~ 0.02 GeV^2
    }
}

TEST_CASE("CCCOH does not apply to hydrogen; cross section rises with E", "[coherent]") {
  YAML::Node c = CohConfig("1000010010:1", 6.0, 14, "none");
  CHECK_THROWS(Generator(c));
  auto sigma = [](double e) {
    Generator gen(CohConfig("1000060120:1", e, 14, "none"));
    return gen.Selector().Get(0).table.sigma[0];
  };
  const double s1 = sigma(1.0), s6 = sigma(6.0), s20 = sigma(20.0);
  CHECK(s1 > 0);
  CHECK(s6 > s1);
  CHECK(s20 > s6);
  // order of magnitude (BS09 Fig. 3: ~0.5e-39 cm^2 / 12C at 1 GeV; RS-type models few e-39 at 6 GeV)
  CHECK(s1 > 0.01);
  CHECK(s1 < 0.15);
  CHECK(s6 > 0.1);
  CHECK(s6 < 0.8);
}
