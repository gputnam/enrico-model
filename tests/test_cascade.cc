// Validation of the intranuclear cascade (theory/fsi_cascade.md):
//  1. 12C density: normalised to A, reproduces the point rms radius, central density ~0.18 fm^-3;
//  2. free cross sections: Cugnon NN values at reference momenta, Delta peak of pi+ p ~ 190-200 mb;
//  3. events: exact conservation of four-momentum, charge and baryon number; every final hadron
//     is physical (status 1, on shell); hydrogen events are untouched;
//  4. hadron-nucleus: pi+ 12C absorption at 165 MeV and proton transparency in the experimental range.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "core/Constants.h"
#include "core/PDG.h"
#include "fsi/HadronNucleonXSec.h"
#include "fsi/IntranuclearCascade.h"
#include "orchestration/Generator.h"
#include "orchestration/OutputWriter.h"

using namespace nuflow;
using Catch::Approx;

namespace {
const Nucleus& Carbon() { return NuclearDatabase::Instance().Get(1000060120); }

struct MemoryWriter : OutputWriter {
  std::vector<Event> events;
  void Begin(const RunInfo&) override {}
  void Write(const Event& e) override { events.push_back(e); }
  void End(const RunSummary&) override {}
};
}  // namespace

TEST_CASE("12C nucleon density", "[cascade]") {
  const NuclearDensity d(Carbon());
  CHECK(d.HarmonicOscillator());
  double sum = 0, r2 = 0;
  const double dr = 0.001;
  for (double r = dr / 2; r < 10; r += dr) {
    sum += 4 * constants::kPi * r * r * d.Rho(r) * dr;
    r2 += 4 * constants::kPi * r * r * r * r * d.Rho(r) * dr;
  }
  CHECK(sum == Approx(12).epsilon(1e-3));
  CHECK(std::sqrt(r2 / sum) == Approx(d.PointRms()).epsilon(2e-3));
  CHECK(d.PointRms() == Approx(2.34).margin(0.02));
  CHECK(d.Rho(0) == Approx(0.179).margin(0.005));
  CHECK(d.KFermi(pdg::kProton, 0) == Approx(0.272).margin(0.005));
  // sampled radii follow r^2 rho
  Random rng(5);
  double s2 = 0;
  const int N = 200000;
  for (int i = 0; i < N; ++i) s2 += d.SamplePosition(rng).Mag2();
  CHECK(std::sqrt(s2 / N) == Approx(d.PointRms()).epsilon(0.01));
}

TEST_CASE("Free hadron-nucleon cross sections", "[cascade]") {
  // Cugnon et al. (1996) at reference momenta; continuity across the piecewise boundaries
  CHECK(hadxs::NNTotal(true, 0.4) == Approx(34.0));
  CHECK(hadxs::NNElastic(true, 2.0) == Approx(22.0).margin(0.2));
  CHECK(hadxs::NNTotal(false, 1.0) == Approx(hadxs::NNTotal(false, 0.9999)).margin(0.1));
  CHECK(hadxs::NNTotal(true, 1.5) == Approx(hadxs::NNTotal(true, 1.4999)).margin(0.1));
  for (double p = 0.1; p < 6; p += 0.01) {
    CHECK(hadxs::NNElastic(true, p) <= hadxs::NNTotal(true, p) + 1e-9);
    CHECK(hadxs::NNElastic(false, p) <= hadxs::NNTotal(false, p) + 1e-9);
  }
  // Delta peak: sigma(pi+ p) = 8 pi / q^2 at W = M_Delta (~195 mb), width 0.11-0.12 GeV
  const double W = hadxs::kMassDelta;
  CHECK(hadxs::DeltaWidth(W) == Approx(0.115).margin(0.005));
  const double q = hadxs::CMMomentum(W, constants::kMassNucleon, constants::kMassPionAverage);
  CHECK(8 * constants::kPi / (q * q) * 0.389379 == Approx(195).margin(10));
  // Oset-Salcedo at T = 165 MeV, rho0: C_Q ~ 13, C_A2 ~ 16, C_A3 ~ 15 MeV
  const auto se = hadxs::DeltaSelfEnergy(0.165, 1.0);
  CHECK(se.cQ * 1e3 == Approx(12.9).margin(0.3));
  CHECK(se.cA2 * 1e3 == Approx(16.3).margin(0.3));
  CHECK(se.cA3 * 1e3 == Approx(15.4).margin(0.4));
}

TEST_CASE("Cascade events conserve energy, momentum, charge and baryon number", "[cascade]") {
  YAML::Node c = YAML::Load(R"(
run: {seed: 7, batch_size: 500, threads: 4}
beam: {pdg: 14, flux: {type: mono, energy: 3.0}}
target: "1000060120:1,1000010010:1"
channels:
  - type: CCQE
    nuclear_model: {type: RFG}
    form_factors: {vector: BBBA05, axial: dipole, MA: 1.03}
  - type: CCRES
    nuclear_model: {type: RFG}
    delta: {C5A0: 1.00, MA: 0.93}
  - type: CCDIS
    nuclear_model: {type: RFG}
    w_min: 1.6
fsi: {type: cascade}
xsec_table: {nodes: 1, samples: 20000, cache: none}
)");
  Generator gen(c);
  MemoryWriter w;
  gen.Run(6000, {&w});
  int nRemnant = 0, nRescattered = 0;
  for (const auto& ev : w.events) {
    FourVector in, out;
    int q = 0, b = 0, qIn = 0, bIn = 0;
    bool hasRemnant = false;
    int nFinalHadrons = 0, nPre = 0;
    for (const auto& p : ev.particles) {
      if (p.status == kStatusBeam) in = in + p.p;
      if (p.status == kStatusTarget) {
        in = in + p.p;
        qIn = pdg::IonZ(p.pdg);
        bIn = pdg::IonA(p.pdg);
      }
      if (p.status == kStatusPreFSI) ++nPre;
      if (p.status == kStatusFinal) {
        out = out + p.p;
        // on-shell final particles
        const double m = p.p.M();
        if (p.pdg == pdg::kProton || p.pdg == pdg::kNeutron) {
          CHECK(m == Approx(pdg::NucleonMass(p.pdg)).margin(1e-6));
          ++b;
          q += p.pdg == pdg::kProton;
          ++nFinalHadrons;
        } else if (std::abs(p.pdg) == 211 || p.pdg == 111) {
          CHECK(m == Approx(pdg::PionMass(p.pdg)).margin(1e-6));
          q += p.pdg == 211 ? 1 : (p.pdg == -211 ? -1 : 0);
          ++nFinalHadrons;
        } else if (p.pdg == 13) {
          q -= 1;
        }
      }
      if (p.status == kStatusRemnantPostFSI) {
        out = out + p.p;
        q += pdg::IonZ(p.pdg);
        b += pdg::IonA(p.pdg);
        hasRemnant = true;
      }
    }
    REQUIRE(std::abs(in.e - out.e) < 1e-8);
    REQUIRE((in.Vect() - out.Vect()).Mag() < 1e-8);
    REQUIRE(q == qIn);
    REQUIRE(b == bIn);
    if (ev.targetPdg == pdg::kProton) {
      CHECK(!hasRemnant);
      CHECK(nFinalHadrons == nPre);  // free nucleon: transparent
    } else {
      nRemnant += hasRemnant;
      nRescattered += nFinalHadrons != nPre;
    }
  }
  CHECK(nRemnant > 0);
  CHECK(nRescattered > 0);
}

TEST_CASE("pi+ 12C absorption and proton transparency", "[cascade]") {
  const Nucleus& C = Carbon();
  IntranuclearCascade inc;
  const auto& d = IntranuclearCascade::DensityFor(C);
  Random rng(11);
  const double R = d.RMax();
  // pi+ at 165 MeV. Data (Ashery et al. PRC 23, 2173 (1981)): sigma_abs = 194 +- 36 mb, sigma_abs+cx = 240 +- 28 mb.
  // The untuned Oset-Salcedo cascade gives about half the absorption at the Delta peak (documented in
  // theory/fsi_cascade.md, Sec. 5); this is a regression check of the model value, not a data comparison.
  const int N = 40000;
  int nAbs = 0, nAbsCx = 0;
  for (int i = 0; i < N; ++i) {
    double bx, by;
    do {
      bx = R * (2 * rng.Uniform() - 1);
      by = R * (2 * rng.Uniform() - 1);
    } while (bx * bx + by * by > R * R);
    IntranuclearCascade::Hadron h;
    h.pdg = pdg::kPiPlus;
    h.x = {bx, by, -R - 0.01};
    h.dir = {0, 0, 1};
    h.estar = 0.165;
    std::vector<IntranuclearCascade::Hadron> esc;
    int z = 6, n = 6;
    inc.Cascade(C, d, {h}, esc, z, n, rng);
    int nPip = 0, nPi = 0;
    for (const auto& e : esc) {
      nPip += e.pdg == 211;
      nPi += std::abs(e.pdg) == 211 || e.pdg == 111;
    }
    nAbs += nPi == 0;
    nAbsCx += nPip == 0;
  }
  const double area = constants::kPi * R * R * 10;
  const double sAbs = area * nAbs / N, sAbsCx = area * nAbsCx / N;
  INFO("sigma_abs = " << sAbs << " mb, sigma_abs+cx = " << sAbsCx << " mb");
  CHECK(sAbs > 80);
  CHECK(sAbs < 110);
  CHECK(sAbsCx > 125);
  CHECK(sAbsCx < 165);

  // Proton transparency at T_p ~ 1-2 GeV: (e,e'p) data ~0.55-0.6 (Garrow et al. PRC 66 044613)
  int clean = 0;
  for (int i = 0; i < N; ++i) {
    IntranuclearCascade::Hadron h;
    h.pdg = pdg::kProton;
    h.x = d.SamplePosition(rng);
    const double c = 2 * rng.Uniform() - 1, s = std::sqrt(1 - c * c), phi = 2 * constants::kPi * rng.Uniform();
    h.dir = {s * std::cos(phi), s * std::sin(phi), c};
    h.estar = 1.5;
    std::vector<IntranuclearCascade::Hadron> esc;
    int z = 5, n = 6;
    IntranuclearCascade::Stats st;
    inc.Cascade(C, d, {h}, esc, z, n, rng, &st);
    clean += st.nucleonCollisions == st.nucleonPauliBlocked;
  }
  const double T = double(clean) / N;
  INFO("transparency(1.5 GeV) = " << T);
  CHECK(T > 0.45);
  CHECK(T < 0.75);
}
