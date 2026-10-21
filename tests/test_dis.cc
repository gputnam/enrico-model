// Tests of the CC inclusive inelastic (DIS) channel: GRV98 LO PDFs, Bodek-Yang structure functions,
// the lepton-hadron tensor contraction, hadronization and end-to-end event generation.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "core/Constants.h"
#include "core/PDG.h"
#include "kinematics/BodekYang.h"
#include "kinematics/DISHadronization.h"
#include "kinematics/DISKinematics.h"
#include "kinematics/GRV98LO.h"
#include "orchestration/Generator.h"

using namespace nuflow;
using namespace nuflow::constants;
using Catch::Approx;

namespace {

// integral over x in [1e-6, 1) on a log grid of f(x) dx
template <typename F>
double IntegrateX(F f) {
  const int n = 20000;
  const double l0 = std::log(1e-6), l1 = 0;
  double s = 0;
  for (int i = 0; i < n; ++i) {
    const double l = l0 + (l1 - l0) * (i + 0.5) / n;
    const double x = std::exp(l);
    if (x >= 1) continue;
    s += f(x) * x * (l1 - l0) / n;
  }
  return s;
}

YAML::Node DisConfig(const std::string& target, double E, int nu = 14, bool withOthers = false) {
  YAML::Node c = YAML::Load(R"(
run: {seed: 7, batch_size: 500, threads: 4}
channels:
  - type: CCDIS
    nuclear_model: {type: RFG}
    w_min: 1.6
fsi: {type: none}
xsec_table: {samples: 100000, cache: none}
)");
  if (withOthers) {
    c["channels"].push_back(YAML::Load("{type: CCQE, nuclear_model: {type: RFG}, form_factors: {vector: BBBA05, axial: dipole, MA: 1.03}}"));
    c["channels"].push_back(YAML::Load("{type: CCRES, nuclear_model: {type: RFG}, delta: {C5A0: 1.0, MA: 0.93}, background: true, nucleon_form_factors: {vector: dipole, axial: dipole, gA: 1.26, MA: 1.05}, w_max: 1.6}"));
  }
  c["beam"]["pdg"] = nu;
  c["beam"]["flux"] = YAML::Load("{type: mono, energy: " + std::to_string(E) + "}");
  c["target"] = target;
  return c;
}

struct MemWriter : OutputWriter {
  std::vector<Event> events;
  void Begin(const RunInfo&) override {}
  void Write(const Event& e) override { events.push_back(e); }
  void End(const RunSummary&) override {}
};

}  // namespace

TEST_CASE("GRV98 LO: momentum and valence sum rules", "[dis]") {
  const auto& g = GRV98LO::Instance();
  for (double Q2 : {2.0, 10.0, 100.0}) {
    const double mom = IntegrateX([&](double x) {
      const auto p = g.Evaluate(x, Q2);
      return p.uv + p.dv + 2 * p.ubar + 2 * p.dbar + 2 * p.s + p.g;
    });
    const double nu = IntegrateX([&](double x) { return g.Evaluate(x, Q2).uv / x; });
    const double nd = IntegrateX([&](double x) { return g.Evaluate(x, Q2).dv / x; });
    INFO("Q2 = " << Q2 << ": momentum " << mom << ", n_u " << nu << ", n_d " << nd);
    CHECK(mom == Approx(1.0).epsilon(0.02));
    CHECK(nu == Approx(2.0).epsilon(0.03));
    CHECK(nd == Approx(1.0).epsilon(0.03));
  }
  // A familiar value: x u_v(0.1, 10 GeV^2) ~ 0.5-0.6 in LO fits.
  const auto p = g.Evaluate(0.1, 10);
  CHECK(p.uv > 0.45);
  CHECK(p.uv < 0.7);
}

TEST_CASE("Tensor contraction reproduces the rest-frame d2sigma/dxdy (with lepton mass)", "[dis]") {
  // Kretzer & Reno, PRD 66, 113007 (2002), Eq. (1) (also Paschos & Yu):
  // d2s/dxdy = G^2 M E / pi { (y^2 x + m^2 y/(2EM)) F1 + [(1 - m^2/4E^2) - (1 + Mx/2E) y] F2
  //            +- [x y (1 - y/2) - m^2 y/(4EM)] F3 + m^2 (m^2 + Q^2)/(4 E^2 M^2 x) F4 - m^2/(EM) F5 }
  const double M = kMassNucleon, m = kMassMuon, E = 5.0;
  StructureFunctions F;
  F.F1 = 0.7;
  F.F2 = 1.1;
  F.F3 = 0.9;
  F.F4 = 0.13;
  F.F5 = 0.8;
  for (int h : {+1, -1})
    for (double x : {0.1, 0.4})
      for (double y : {0.2, 0.6}) {
        const double nu = y * E, Q2 = 2 * M * nu * x;
        const double Ep = E - nu, pp = std::sqrt(Ep * Ep - m * m);
        const double cosT = (2 * E * Ep - m * m - Q2) / (2 * E * pp);
        const FourVector k(E, 0, 0, E), kp(Ep, pp * std::sqrt(1 - cosT * cosT), 0, pp * cosT), P(M, 0, 0, 0);
        const FourVector q = k - kp;
        const double mine = kGF * kGF * M * nu / (16 * kPi * E) * DISKinematics::Contract(k, kp, P, M, q, F, h);
        const double ref = kGF * kGF * M * E / kPi *
                           ((y * y * x + m * m * y / (2 * E * M)) * F.F1 +
                            ((1 - m * m / (4 * E * E)) - (1 + M * x / (2 * E)) * y) * F.F2 +
                            h * (x * y * (1 - y / 2) - m * m * y / (4 * E * M)) * F.F3 +
                            m * m * (m * m + Q2) / (4 * E * E * M * M * x) * F.F4 - m * m / (E * M) * F.F5);
        INFO("h=" << h << " x=" << x << " y=" << y);
        CHECK(mine == Approx(ref).epsilon(1e-9));
      }
}

TEST_CASE("Bodek-Yang structure functions: limits", "[dis]") {
  BodekYang by;
  const double M = kMassNucleon;
  // High Q2: K factors -> 1, xi_w -> x; F2(nu p) ~ 2x(d + s + ubar), Callan-Gross 2xF1 ~ F2/(1+R).
  const double x = 0.3, Q2 = 200;
  const auto F = by.Compute(true, true, x, Q2, Q2 / (2 * M * x), M);
  const double xi = by.XiW(x, Q2, M, 0);
  CHECK(xi == Approx(x).epsilon(0.01));
  const auto p = GRV98LO::Instance().Evaluate(xi, Q2);
  CHECK(F.F2 == Approx(2 * 1.026 * (p.dv + p.dbar + p.ubar + p.s)).epsilon(0.03));
  CHECK(2 * x * F.F1 / F.F2 == Approx(1 / (1 + BodekYang::R(x, Q2))).epsilon(0.02));
  CHECK(F.F3 > 0);
  CHECK(F.F5 == Approx(F.F2 / (2 * x)));
  // Q2 -> 0 at fixed nu: the vector part vanishes, the axial (PCAC) part of F2 stays finite, so
  // 2xF1 (= 2 x 2xF1^V) / F2 -> 0.
  const double nu0 = 2.0, q20 = 1e-4, x0 = q20 / (2 * M * nu0);
  const auto F0 = by.Compute(true, true, x0, q20, nu0, M);
  CHECK(F0.F2 > 0);
  CHECK(2 * x0 * F0.F1 < 1e-2 * F0.F2);
  // R1998 values of order 0.2-0.3 at x ~ 0.1, Q2 ~ 2 (E143 Fig. 3/4).
  CHECK(BodekYang::R1998(0.1, 2.0) > 0.15);
  CHECK(BodekYang::R1998(0.1, 2.0) < 0.4);
}

TEST_CASE("Hadronization: conservation laws and charged multiplicity", "[dis]") {
  DISHadronization had;
  Random rng(11);
  for (int nuPdg : {14, -14})
    for (int hit : {pdg::kProton, pdg::kNeutron})
      for (double W : {1.7, 2.2, 3.0, 5.0, 8.0}) {
        const FourVector X(std::sqrt(W * W + 9.0), 1.0, 0.5, 2.7);  // arbitrary lab momentum
        const int QX = (hit == pdg::kProton ? 1 : 0) + (nuPdg > 0 ? 1 : -1);
        double sn = 0, sn2 = 0;
        const int N = 4000;
        for (int i = 0; i < N; ++i) {
          std::vector<DISHadronization::Hadron> out;
          had.Hadronize(nuPdg, hit, X, {0.3, 0.1, 1.0}, rng, out);
          REQUIRE(!out.empty());
          FourVector sum;
          int Q = 0, B = 0, nch = 0, nMes = 0;
          for (const auto& h : out) {
            sum = sum + h.p;
            const int c = h.pdg == pdg::kProton || h.pdg == pdg::kPiPlus ? 1 : (h.pdg == pdg::kPiMinus ? -1 : 0);
            Q += c;
            nch += c != 0;
            B += (h.pdg == pdg::kProton || h.pdg == pdg::kNeutron);
            nMes += (h.pdg == pdg::kPiPlus || h.pdg == pdg::kPiMinus || h.pdg == pdg::kPiZero);
            const double m = (h.pdg == pdg::kProton || h.pdg == pdg::kNeutron) ? pdg::NucleonMass(h.pdg) : pdg::PionMass(h.pdg);
            REQUIRE(std::abs(h.p.M() - m) < 1e-6);
          }
          REQUIRE(Q == QX);
          REQUIRE(B == 1);
          REQUIRE(nMes >= 1);
          REQUIRE(std::abs(sum.e - X.e) < 1e-7);
          REQUIRE((sum.Vect() - X.Vect()).Mag() < 1e-7);
          sn += nch;
          sn2 += nch * nch;
        }
        const double mean = sn / N, D = std::sqrt(sn2 / N - mean * mean);
        const double expected = DISHadronization::MeanChargedMultiplicity(nuPdg, hit, W);
        INFO("nu " << nuPdg << " hit " << hit << " W " << W << ": <nch> " << mean << " (data fit " << expected
                   << "), D/<nch> " << D / mean);
        if (W >= 3.0) {
          CHECK(mean == Approx(expected).epsilon(0.08));
          CHECK(D / mean == Approx(0.35).margin(0.07));
        }
      }
}

TEST_CASE("Free-nucleon DIS cross section vs world isoscalar sigma/E", "[dis]") {
  // Isoscalar nucleon (p + n)/2 at 40 GeV: the world average sigma/E = 0.675e-38 cm^2/GeV (30-50 GeV,
  // PDG; quoted by Bodek-Yang). The W > 1.6 GeV part carries ~95% of it at this energy.
  // Free proton, and carbon per nucleon as the isoscalar proxy (binding and Fermi motion change the
  // W > 1.6 GeV integral by ~1.5% at 40 GeV, theory/dis.md Sec. 1.3).
  Generator gp(DisConfig("1000010010:1", 40.0));
  const double sp = gp.Selector().Get(0).table.sigma[0];
  Generator gc(DisConfig("1000060120:1", 40.0));
  const double sc = gc.Selector().Get(0).table.sigma[0] / 12;
  INFO("sigma/E: proton " << sp / 40 << ", carbon per nucleon " << sc / 40);
  CHECK(sc / 40 > 0.60);
  CHECK(sc / 40 < 0.72);
  CHECK(sp / sc > 0.5);  // nu p < nu n
  CHECK(sp / sc < 0.8);
}

TEST_CASE("CCDIS events: four-momentum conservation and final-state content", "[dis]") {
  for (int nu : {14, -14}) {
    Generator gen(DisConfig("1000060120:1,1000010010:1", 6.0, nu));
    MemWriter w;
    gen.Run(3000, {&w});
    for (const auto& ev : w.events) {
      FourVector in, out;
      int nLep = 0, nBar = 0, nMes = 0, Q = 0, Qin = 0;
      for (const auto& p : ev.particles) {
        if (p.status == kStatusBeam || p.status == kStatusTarget) in = in + p.p;
        if (p.status == kStatusFinal || p.status == kStatusRemnant) out = out + p.p;
        if (p.status == kStatusStruckNucleon) Qin = p.pdg == pdg::kProton ? 1 : 0;
        if (p.status != kStatusFinal) continue;
        const int a = std::abs(p.pdg);
        if (a == 13) {
          ++nLep;
          Q += p.pdg == 13 ? -1 : 1;
        } else if (a == 2212 || a == 2112) {
          ++nBar;
          Q += p.pdg == 2212;
        } else if (a == 211 || a == 111) {
          ++nMes;
          Q += p.pdg == 211 ? 1 : (p.pdg == -211 ? -1 : 0);
        }
      }
      CHECK(nLep == 1);
      CHECK(nBar == 1);
      CHECK(nMes >= 1);
      CHECK(Q == Qin);
      CHECK(std::abs(in.e - out.e) < 1e-7);
      CHECK((in.Vect() - out.Vect()).Mag() < 1e-7);
      CHECK(ev.W >= 1.6 - 1e-12);
      CHECK(ev.channel == 600);
    }
  }
}
