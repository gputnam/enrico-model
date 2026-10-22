// Validation of the CCNSTAR channel (N(1440), N(1520), N(1535); theory/nstar.md):
//  1. the numerical Dirac-trace hadron tensor reproduces the closed-form structure functions of
//     Lalakulich, Paschos, Piranishvili, PRD 74, 014009 (2006), Eqs. (4.6)-(4.8) (spin 3/2, lower
//     sign) and (4.11)-(4.13) (spin 1/2), for nucleon at rest, massless lepton, nu and nubar;
//  2. partial widths and Breit-Wigner normalisation;
//  3. events conserve four-momentum and charge, have one lepton, one nucleon and >= 1 meson, W < w_max,
//     and only nu n / nubar p are excited.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "core/Constants.h"
#include "core/PDG.h"
#include "kinematics/NStarKinematics.h"
#include "orchestration/Generator.h"

using namespace nuflow;
using Catch::Approx;

namespace {

// LPP structure functions V1, V2, V3 (with C3A = C4A = 0 for spin 3/2), W_i = c V_i.
void LPPV(const NStarResonance& r, const NStarFormFactors& f, double Q2, double pq, double V[3]) {
  const double m = constants::kMassNucleon, M = r.M, m2 = m * m, M2 = M * M;
  if (r.twoJ == 3) {
    const double s = r.parity > 0 ? 1.0 : -1.0;  // upper sign: P33, lower: D13
    const double C3V = f.C3V, C4V = f.C4V, C5V = f.C5V, C5A = f.C5A;
    const double a = Q2 - pq;
    V[0] = C3V * C3V / (m2 * M2) * (a * a * (pq + m2) + M2 * (pq * pq + Q2 * m2 + s * Q2 * m * M)) +
           (C3V * C4V * a - C3V * C5V * pq) / (m2 * m * M) * (a * (pq + m2 - s * 2 * m * M) - M2 * pq) +
           std::pow(C4V * a - C5V * pq, 2) / (m2 * m2) * (pq + m2 - s * m * M) + C5A * C5A * (pq + m2 + s * m * M);
    V[1] = C3V * C3V / M2 * Q2 * (pq + m2 + M2) + C3V * C4V / (m * M) * Q2 * (pq + std::pow(M - s * m, 2)) +
           C3V * C5V / (m * M) * Q2 * (pq + std::pow(M - s * m, 2) + Q2) +
           (C4V * C4V / m2 + C5V * C5V * (Q2 + M2) / (m2 * M2) + 2 * C4V * C5V / m2) * Q2 * (pq + m2 - s * m * M) +
           C5A * C5A * m2 / M2 * (pq + m2 + s * m * M);
    V[2] = -2 * C4V * C5A * a + 2 * C5V * C5A * pq + 2 * C3V * C5A * m / M * (2 * M2 + s * 2 * m * M + Q2 - pq);
    for (int i = 0; i < 3; ++i) V[i] *= 2.0 / (3 * m);
  } else {
    const double s = r.parity > 0 ? 1.0 : -1.0;  // upper sign: P11, lower: S11
    const double mu = m + M, g1V = f.g1V, g2V = f.g2V, g1A = f.g1A;
    V[0] = g1V * g1V / std::pow(mu, 4) * Q2 * Q2 * (pq + m2 - s * m * M) +
           g2V * g2V / (mu * mu) * (2 * pq * pq + Q2 * (m2 + s * m * M - pq)) +
           g1V * g2V / std::pow(mu, 3) * 2 * Q2 * (pq * (M - s * m) + s * m * Q2) + g1A * g1A * (m2 + s * m * M + pq);
    V[1] = 2 * m2 * (g1V * g1V / std::pow(mu, 4) * Q2 * Q2 + g2V * g2V / (mu * mu) * Q2 + g1A * g1A);
    V[2] = 4 * m2 * (g1V * g1A / (mu * mu) * Q2 + g2V * g1A / mu * (M + s * m));
    for (int i = 0; i < 3; ++i) V[i] *= 1.0 / m;
  }
}

struct MemWriter : OutputWriter {
  std::vector<Event> events;
  void Begin(const RunInfo&) override {}
  void Write(const Event& e) override { events.push_back(e); }
  void End(const RunSummary&) override {}
};

YAML::Node NStarConfig(const std::string& target, double e, int nu) {
  YAML::Node c = YAML::Load(R"(
run: {seed: 7, batch_size: 500, threads: 4}
channels:
  - type: CCNSTAR
    nuclear_model: {type: RFG}
fsi: {type: none}
xsec_table: {nodes: 2, samples: 100000, cache: none}
)");
  c["beam"]["pdg"] = nu;
  c["beam"]["flux"]["type"] = "mono";
  c["beam"]["flux"]["energy"] = e;
  c["target"] = target;
  return c;
}

}  // namespace

TEST_CASE("N* hadron tensor reproduces the LPP closed-form structure functions", "[nstar]") {
  NStarKinematics kin(NStarKinematics::Options{});
  const double m = constants::kMassNucleon;
  const FourVector p(m, 0, 0, 0);
  for (const auto& r : kin.Resonances()) {
    for (double E : {1.5, 3.0, 8.0})
      for (double Q2 : {0.05, 0.4, 1.2}) {
        const double W = r.M;
        const double q0 = (W * W - m * m + Q2) / (2 * m);
        const double Ep = E - q0;
        if (Ep <= 0) continue;
        const double cosT = 1 - Q2 / (2 * E * Ep);
        if (cosT < -1) continue;
        const double sinT = std::sqrt(1 - cosT * cosT);
        const FourVector k(E, 0, 0, E), kp(Ep, Ep * sinT, 0, Ep * cosT);
        const FourVector pR = p + k - kp;
        REQUIRE(pR.M() == Approx(W).epsilon(1e-9));
        const double pq = m * q0, kpd = m * E, kpp = m * Ep;
        double V[3];
        LPPV(r, kin.FormFactors(r, Q2), Q2, pq, V);
        for (int hel : {+1, -1}) {
          const double lpp = V[0] * Q2 + V[1] / (m * m) * (2 * kpd * kpp - 0.5 * m * m * Q2) + hel * V[2] / (m * m) * (Q2 * kpd - 0.5 * pq * Q2);
          const double lh = kin.ContractLH(r, k, kp, 0.0, hel, p, m, pR);
          INFO(r.name << " E=" << E << " Q2=" << Q2 << " hel=" << hel);
          CHECK(lh == Approx(32 * m * lpp).epsilon(1e-8));
        }
      }
  }
}

TEST_CASE("N* widths and Breit-Wigner", "[nstar]") {
  NStarKinematics kin(NStarKinematics::Options{});
  for (const auto& r : kin.Resonances()) {
    double a, b, c;
    kin.PartialWidths(r, r.M, a, b, c);
    CHECK(a + b + c == Approx(r.Gamma0).epsilon(1e-9));
    CHECK(a == Approx(r.xPi * r.Gamma0).epsilon(1e-9));
    kin.PartialWidths(r, 1.10, a, b, c);
    CHECK(b == 0.0);
    CHECK(c == 0.0);  // below the N pi pi threshold
    // integral of A(W^2) dW^2 up to W = 2 GeV: 0.6-1 (the P11 is broad; the energy-dependent widths
    // grow above the pole, so the far tail of the D-wave D13 is enhanced)
    double sum = 0;
    for (double W = 1.0785; W < 2.0; W += 0.0005) sum += kin.Spectral(r, W) * 2 * W * 0.0005;
    CHECK(sum > 0.6);
    CHECK(sum < 1.0);
  }
}

TEST_CASE("CCNSTAR events: conservation laws, final-state content, isospin selection", "[nstar]") {
  for (int nu : {14, -14}) {
    // nu on hydrogen has no I = 1/2 CC channel: carbon only for neutrinos.
    Generator gen(NStarConfig(nu > 0 ? "1000060120:1" : "1000060120:1,1000010010:1", 3.0, nu));
    MemWriter w;
    gen.Run(2000, {&w});
    int nEta = 0, nMulti = 0;
    for (const auto& ev : w.events) {
      FourVector in, out;
      int nLep = 0, nNuc = 0, nMes = 0, charge = 0, hit = 0;
      for (const auto& p : ev.particles) {
        if (p.status == kStatusBeam || p.status == kStatusTarget) in = in + p.p;
        if (p.status == kStatusFinal || p.status == kStatusRemnant) out = out + p.p;
        if (p.status == kStatusStruckNucleon) hit = p.pdg;
        if (p.status != kStatusFinal) continue;
        const int a = std::abs(p.pdg);
        if (a == 13) ++nLep, charge += p.pdg > 0 ? -1 : 1;
        else if (a == 2212 || a == 2112) ++nNuc, charge += p.pdg == 2212;
        else if (a == 211 || a == 111 || a == 221) ++nMes, charge += p.pdg == 211 ? 1 : p.pdg == -211 ? -1 : 0;
        if (a == 221) ++nEta;
      }
      if (nMes >= 2) ++nMulti;
      CHECK(hit == (nu > 0 ? 2112 : 2212));
      CHECK(charge == (nu > 0 ? 0 : 1));  // hadrons + lepton = struck nucleon charge
      CHECK(nLep == 1);
      CHECK(nNuc == 1);
      CHECK(nMes >= 1);
      CHECK(std::abs(in.e - out.e) < 1e-9);
      CHECK(std::abs(in.pz - out.pz) < 1e-9);
      CHECK(std::abs(in.px - out.px) < 1e-9);
      CHECK(ev.W <= 1.6 + 1e-12);
      CHECK(ev.channel == 401);
    }
    CHECK(nEta > 0);
    CHECK(nMulti > 0);
  }
}

namespace {
YAML::Node DisIsoConfig(const std::string& target, double e, bool iso) {
  YAML::Node c = YAML::Load(R"(
run: {seed: 9, batch_size: 500, threads: 4}
channels:
  - type: CCDIS
    nuclear_model: {type: RFG}
    w_min: 1.6
    w_max: 2.0
fsi: {type: none}
xsec_table: {nodes: 2, samples: 100000, cache: none}
)");
  if (iso) c["channels"][0]["single_pion"] = YAML::Load("{scale_i32: 0.28, scale_i12: 0.68, isospin: true}");
  c["beam"]["pdg"] = 14;
  c["beam"]["flux"]["type"] = "mono";
  c["beam"]["flux"]["energy"] = e;
  c["target"] = target;
  return c;
}
}  // namespace

TEST_CASE("CCDIS single-pion isospin decomposition: inclusive unchanged, 1pi share scaled", "[nstar][dis]") {
  for (const char* tgt : {"1000010010:1", "1000000010:1"}) {
    Generator g0(DisIsoConfig(tgt, 4.0, false)), g1(DisIsoConfig(tgt, 4.0, true));
    const double s0 = g0.Selector().Get(0).table.sigma[0], s1 = g1.Selector().Get(0).table.sigma[0];
    CHECK(s1 == Approx(s0).epsilon(0.01));  // weights unchanged (composition only); tables differ by MC noise
    double f[2];
    Generator* gs[2] = {&g0, &g1};
    for (int i = 0; i < 2; ++i) {
      MemWriter w;
      gs[i]->Run(4000, {&w});
      int single = 0;
      for (const auto& ev : w.events) {
        int nPi = 0, charge = 0;
        FourVector in, out;
        for (const auto& p : ev.particles) {
          if (p.status == kStatusBeam || p.status == kStatusTarget) in = in + p.p;
          if (p.status == kStatusFinal || p.status == kStatusRemnant) out = out + p.p;
          if (p.status != kStatusFinal) continue;
          if (p.pdg == 211 || p.pdg == -211 || p.pdg == 111) ++nPi;
          charge += p.pdg == 2212 || p.pdg == 211 ? 1 : (p.pdg == 13 || p.pdg == -211) ? -1 : 0;
        }
        CHECK(std::abs(in.e - out.e) < 1e-9);
        CHECK(charge == (std::string(tgt) == "1000010010:1" ? 1 : 0));
        CHECK(nPi >= 1);
        single += nPi == 1;
      }
      f[i] = double(single) / w.events.size();
    }
    if (std::string(tgt) == "1000010010:1") CHECK(f[1] == Approx(0.28 * f[0]).margin(0.03));  // pure I = 3/2
    else CHECK(f[1] < f[0]);
  }
}
