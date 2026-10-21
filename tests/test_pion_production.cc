// Validation of the HNV pion-production current and the CCRES channel:
//  1. the optimised Delta-pole amplitudes == a literal transcription of HNV Eqs. (40), (44), (51);
//  2. CVC: q_mu V^mu = 0 between on-shell spinors, for every isospin channel, with background;
//  3. Contract() reproduces the closed-form CCQE contraction for a nucleon current (conventions, nu/nubar);
//  4. Delta-only isospin ratios 1 : 1/9 : 2/9 (p pi+ : n pi+ : p pi0);
//  5. Delta width at the pole;
//  6. free-proton total cross section vs ANL/BNL bubble-chamber data; events conserve four-momentum
//     and have exactly one lepton, one nucleon and one pion in the final state.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "core/Constants.h"
#include "core/Dirac.h"
#include "core/PDG.h"
#include "core/Random.h"
#include "kinematics/PionProductionCurrent.h"
#include "kinematics/QECurrent.h"
#include "orchestration/Generator.h"

using namespace nuflow;
using namespace nuflow::dirac;
using Catch::Approx;

namespace {

FourVector OnShell(double M, const ThreeVector& p) { return {std::sqrt(p.Mag2() + M * M), p}; }
ThreeVector RandomVec(Random& r, double scale) {
  return {scale * (2 * r.Uniform() - 1), scale * (2 * r.Uniform() - 1), scale * (2 * r.Uniform() - 1)};
}
double Lower(const FourVector& v, int mu) { return Metric(mu) * Comp(v, mu); }

NucleonFormFactors HNVNucleonFF() {
  return NucleonFormFactors::FromConfig(YAML::Load("{vector: dipole, axial: dipole, gA: 1.26, MA: 1.05}"));
}

// HNV Eq. (40): Gamma^{alpha mu}(p, q), upper indices, pD = p + q.
Mat GammaRef(int al, int mu, const FourVector& p, const FourVector& q, const NDeltaFormFactors& C) {
  const Gammas& G = Gammas::Get();
  const double M = constants::kMassNucleon, M2 = M * M;
  const FourVector pD = p + q;
  const double g = (al == mu) ? Metric(mu) : 0.0;
  const Mat qs = G.Slash(q);
  Mat V = (qs * g - G.g[mu] * Comp(q, al)) * (C.C3V / M);
  V.AddScaled(G.one, C.C4V / M2 * (g * q.Dot(pD) - Comp(q, al) * Comp(pD, mu)));
  V.AddScaled(G.one, C.C5V / M2 * (g * q.Dot(p) - Comp(q, al) * Comp(p, mu)));
  Mat A = (qs * g - G.g[mu] * Comp(q, al)) * (C.C3A / M);
  A.AddScaled(G.one, C.C4A / M2 * (g * q.Dot(pD) - Comp(q, al) * Comp(pD, mu)));
  A.AddScaled(G.one, C.C5A * g);
  A.AddScaled(G.one, C.C6A / M2 * Comp(q, al) * Comp(q, mu));
  return V * G.g5 + A;
}

// HNV Eq. (44) with lower indices.
Mat ProjRef(int al, int be, const FourVector& pD, double MD) {
  const Gammas& G = Gammas::Get();
  const Mat ga = G.g[al] * Metric(al), gb = G.g[be] * Metric(be);
  Mat B = G.one * ((al == be ? Metric(al) : 0.0) - 2 * Lower(pD, al) * Lower(pD, be) / (3 * MD * MD));
  B.AddScaled(ga * gb, -1.0 / 3);
  B.AddScaled(gb, Lower(pD, al) / (3 * MD));
  B.AddScaled(ga, -Lower(pD, be) / (3 * MD));
  return (G.Slash(pD) + G.one * MD) * B * (-1.0);
}

void ReferenceDelta(const PionProductionCurrent& cur, const FourVector& p, const FourVector& pp, const FourVector& kpi,
                    Mat dp[4], Mat cdp[4]) {
  const Gammas& G = Gammas::Get();
  const cplx I(0, 1);
  const auto& par = cur.Parameters();
  const double mpi = constants::kMassPionAverage, MD = par.MDelta;
  const FourVector q = pp + kpi - p;
  const NDeltaFormFactors C = cur.DeltaFormFactors(-q.M2());
  const FourVector pD = p + q, pDc = pp - q;
  const cplx facD = I * (par.fStar / mpi) * std::sqrt(3.0) / cplx(pD.M2() - MD * MD, MD * cur.DeltaWidth(pD.M2()));
  const cplx facC = I * (par.fStar / mpi) / std::sqrt(3.0) / cplx(pDc.M2() - MD * MD, MD * cur.DeltaWidth(pDc.M2()));
  const FourVector mq(-q.e, -q.px, -q.py, -q.pz);
  for (int mu = 0; mu < 4; ++mu) {
    Mat d, c;
    for (int a = 0; a < 4; ++a)
      for (int b = 0; b < 4; ++b) {
        d += ProjRef(a, b, pD, MD) * GammaRef(b, mu, p, q, C) * Comp(kpi, a);
        c += G.Bar(GammaRef(a, mu, pp, mq, C)) * ProjRef(a, b, pDc, MD) * Comp(kpi, b);
      }
    dp[mu] = d * facD;
    cdp[mu] = c * facC;
  }
}

double MaxAbs(const Mat& m) {
  double r = 0;
  for (const auto& x : m.m) r = std::max(r, std::abs(x));
  return r;
}

// Random on-shell N pi final state with initial nucleon p.
void RandomKinematics(Random& rng, FourVector& p, FourVector& pp, FourVector& kpi) {
  const double M = constants::kMassNucleon, mpi = constants::kMassPionAverage;
  p = OnShell(M, RandomVec(rng, 0.25));
  pp = OnShell(M, RandomVec(rng, 0.8));
  kpi = OnShell(mpi, RandomVec(rng, 0.5));
}

}  // namespace

TEST_CASE("HNV Delta-pole amplitudes match the literal Eqs. (40), (44), (51)", "[res]") {
  PionProductionCurrent cur(HNVParameters{}, HNVNucleonFF());
  Random rng(7);
  for (int trial = 0; trial < 20; ++trial) {
    FourVector p, pp, kpi;
    RandomKinematics(rng, p, pp, kpi);
    PionCurrentPieces P;
    cur.Pieces(p, pp, kpi, P);
    Mat dp[4], cdp[4];
    ReferenceDelta(cur, p, pp, kpi, dp, cdp);
    for (int mu = 0; mu < 4; ++mu) {
      const double scale = std::max(MaxAbs(dp[mu]), 1e-12);
      CHECK(MaxAbs(P.dp[mu] - dp[mu]) / scale < 1e-10);
      const double scaleC = std::max(MaxAbs(cdp[mu]), 1e-12);
      CHECK(MaxAbs(P.cdp[mu] - cdp[mu]) / scaleC < 1e-10);
    }
  }
}

TEST_CASE("Vector current is conserved between on-shell spinors (all isospin channels)", "[res]") {
  HNVParameters par;
  par.axialScale = 0;  // vector current only
  PionProductionCurrent cur(par, HNVNucleonFF());
  const Gammas& G = Gammas::Get();
  const double M = constants::kMassNucleon;
  Random rng(11);
  for (int trial = 0; trial < 20; ++trial) {
    FourVector p, pp, kpi;
    RandomKinematics(rng, p, pp, kpi);
    const FourVector q = pp + kpi - p;
    PionCurrentPieces P;
    cur.Pieces(p, pp, kpi, P);
    for (auto iso : {PionIsospin::PPiPlus, PionIsospin::NPiPlus, PionIsospin::PPiZero}) {
      Mat S[4];
      cur.Combine(P, iso, S);
      Mat div;
      double ref = 0;
      for (int mu = 0; mu < 4; ++mu) {
        const Mat sand = (G.Slash(pp) + G.one * M) * S[mu] * (G.Slash(p) + G.one * M);
        div.AddScaled(sand, Lower(q, mu));
        ref += MaxAbs(sand) * std::abs(Comp(q, mu));
      }
      CHECK(ref > 0);
      CHECK(MaxAbs(div) < 1e-10 * ref);
    }
  }
}

TEST_CASE("Contract() reproduces the closed-form CCQE contraction", "[res]") {
  const Gammas& G = Gammas::Get();
  const cplx I(0, 1);
  Random rng(5);
  const double M = constants::kMassNucleon, m = constants::kMassMuon;
  for (int trial = 0; trial < 20; ++trial) {
    const WeakFormFactors ff{rng.Uniform(-1, 2), rng.Uniform(-1, 4), rng.Uniform(-2, 2), rng.Uniform(-5, 20)};
    const FourVector k = OnShell(0, RandomVec(rng, 3)), kp = OnShell(m, RandomVec(rng, 2));
    const FourVector p = OnShell(M, RandomVec(rng, 0.3)), pp = OnShell(M, RandomVec(rng, 1.5));
    const FourVector q = pp - p;
    Mat S[4];
    for (int mu = 0; mu < 4; ++mu) {
      S[mu] = G.g[mu] * ff.F1;
      S[mu].AddScaled((G.g[mu] * G.Slash(q) - G.Slash(q) * G.g[mu]), -ff.F2 / (4 * M));
      S[mu].AddScaled(G.g5g[mu], -ff.FA);
      S[mu].AddScaled(G.g5, -Comp(q, mu) * ff.FP / M);
    }
    for (int h : {+1, -1})
      CHECK(PionProductionCurrent::Contract(k, kp, m, h, p, M, pp, M, S) ==
            Approx(ContractQE(k, kp, p, pp, ff, M, h)).epsilon(1e-9));
  }
}

TEST_CASE("Delta-only isospin ratios and Delta width", "[res]") {
  HNVParameters par;
  par.background = false;
  PionProductionCurrent cur(par, HNVNucleonFF());
  Random rng(3);
  const double m = constants::kMassMuon, M = constants::kMassNucleon;
  for (int trial = 0; trial < 10; ++trial) {
    FourVector p, pp, kpi;
    RandomKinematics(rng, p, pp, kpi);
    const FourVector k(3.0, 0, 0, 3.0);
    const FourVector kp = k + p - pp - kpi;  // off-shell lepton is fine for a ratio test
    PionCurrentPieces P;
    cur.Pieces(p, pp, kpi, P);
    Mat a[4], b[4], c[4];
    cur.Combine(P, PionIsospin::PPiPlus, a);
    cur.Combine(P, PionIsospin::NPiPlus, b);
    cur.Combine(P, PionIsospin::PPiZero, c);
    const double la = PionProductionCurrent::Contract(k, kp, m, 1, p, M, pp, M, a);
    CHECK(PionProductionCurrent::Contract(k, kp, m, 1, p, M, pp, M, b) / la == Approx(1.0 / 9).epsilon(1e-9));
    CHECK(PionProductionCurrent::Contract(k, kp, m, 1, p, M, pp, M, c) / la == Approx(2.0 / 9).epsilon(1e-9));
  }
  // HNV Eq. (45) at the pole: ~ 0.11-0.12 GeV (PDG Breit-Wigner width 0.114-0.120 GeV).
  const double g = cur.DeltaWidth(1.232 * 1.232);
  CHECK(g > 0.105);
  CHECK(g < 0.125);
  CHECK(cur.DeltaWidth(1.07 * 1.07) == 0.0);
}

namespace {
struct MemWriter : OutputWriter {
  std::vector<Event> events;
  void Begin(const RunInfo&) override {}
  void Write(const Event& e) override { events.push_back(e); }
  void End(const RunSummary&) override {}
};

YAML::Node ResConfig(const std::string& target, double e, int nu = 14) {
  YAML::Node c = YAML::Load(R"(
run: {seed: 5, batch_size: 500, threads: 4}
channels:
  - type: CCRES
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

TEST_CASE("Free-proton nu_mu p -> mu- p pi+ cross section vs ANL/BNL", "[res]") {
  // Wilkinson et al., PRD 90, 112017 (2014), Fig. 6 (ANL+BNL, re-extracted, no W cut):
  // ~0.45-0.6 at 1 GeV, ~0.65-0.8 at 2-3 GeV [1e-38 cm^2].
  Generator g1(ResConfig("1000010010:1", 1.0));
  const double s1 = g1.Selector().Get(0).table.sigma[0];
  CHECK(s1 > 0.35);
  CHECK(s1 < 0.75);
  Generator g3(ResConfig("1000010010:1", 2.5));
  const double s3 = g3.Selector().Get(0).table.sigma[0];
  CHECK(s3 > 0.5);
  CHECK(s3 < 1.0);
}

TEST_CASE("CCRES events: four-momentum conservation and final-state content", "[res]") {
  for (int nu : {14, -14}) {
    Generator gen(ResConfig("1000060120:1,1000010010:1", 3.0, nu));
    MemWriter w;
    gen.Run(2000, {&w});
    for (const auto& ev : w.events) {
      FourVector in, out;
      int nLep = 0, nNuc = 0, nPi = 0;
      for (const auto& p : ev.particles) {
        if (p.status == kStatusBeam || p.status == kStatusTarget) in = in + p.p;
        if (p.status == kStatusFinal || p.status == kStatusRemnant) out = out + p.p;
        if (p.status != kStatusFinal) continue;
        const int a = std::abs(p.pdg);
        if (a == 13) ++nLep;
        else if (a == 2212 || a == 2112) ++nNuc;
        else if (a == 211 || a == 111) ++nPi;
      }
      CHECK(nLep == 1);
      CHECK(nNuc == 1);
      CHECK(nPi == 1);
      CHECK(std::abs(in.e - out.e) < 1e-9);
      CHECK(std::abs(in.pz - out.pz) < 1e-9);
      CHECK(std::abs(in.px - out.px) < 1e-9);
      CHECK(ev.W <= 1.6 + 1e-12);
      CHECK(ev.channel == 400);
    }
  }
}
