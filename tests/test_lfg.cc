// Local Fermi gas and Valencia RPA (theory/local_fermi_gas.md):
//  1. LFG sampling: local Fermi sphere, removal energies, energy gap and outgoing potential for 12C;
//  2. Lindhard function: static and high-frequency limits, imaginary part vs Monte Carlo phase space;
//  3. RPA tensor contraction: per-block factors vs a brute-force Dirac-trace tensor decomposition;
//  4. end-to-end CCQE events with LFG / SF (+RPA) and the cascade: four-momentum conservation, vertex;
//  5. 12C spectral function sampling (theory/spectral_function.md): shells, SRC pairs, energy balance.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "DiracAlgebra.h"
#include "core/Constants.h"
#include "core/Nucleus.h"
#include "core/PDG.h"
#include "core/Random.h"
#include "fsi/NuclearDensity.h"
#include "initial_state/LocalFermiGas.h"
#include "initial_state/ShellModelSpectralFunction.h"
#include "kinematics/QECurrent.h"
#include "kinematics/ValenciaRPA.h"
#include "orchestration/Generator.h"

using namespace nuflow;
using Catch::Approx;

namespace {
constexpr double kHbarC = 0.1973269804;
FourVector OnShell(double M, const ThreeVector& p) { return {std::sqrt(p.Mag2() + M * M), p}; }
ThreeVector RandomVec(Random& r, double scale) {
  return {scale * (2 * r.Uniform() - 1), scale * (2 * r.Uniform() - 1), scale * (2 * r.Uniform() - 1)};
}
}  // namespace

TEST_CASE("Local Fermi gas sampling and energy balance for 12C", "[lfg]") {
  const Nucleus& C = NuclearDatabase::Instance().Get(1000060120);
  const NuclearDensity& d = NuclearDensity::For(C);
  LocalFermiGas lfg;
  Random rng(11);
  HitNucleonRequest req;
  req.species[0] = pdg::kNeutron;
  const int N = 200000;
  double sumER = 0, sumK2 = 0, sumR2 = 0;
  for (int i = 0; i < N; ++i) {
    InitialState s;
    lfg.Sample(C, req, s, rng);
    REQUIRE(s.hasPosition);
    const double r = s.position.Mag();
    const double k = s.hit[0].p.Mag();
    REQUIRE(k <= d.KFermi(pdg::kNeutron, r) + 1e-12);
    // E_R = T_F(r) + S_n - T(k) in [S_n, S_n + T_F(r)]
    REQUIRE(s.hit[0].removalEnergy >= C.separationNeutron - 1e-12);
    const FourVector tot = s.hit[0].Momentum() + s.remnant;
    REQUIRE(std::abs(tot.e - C.mass) < 1e-9);
    sumER += s.hit[0].removalEnergy;
    sumK2 += k * k;
    sumR2 += r * r;
  }
  // The vertex follows the point-nucleon density: <r^2> = rms^2.
  CHECK(std::sqrt(sumR2 / N) == Approx(d.PointRms()).epsilon(0.01));
  // <k^2> = (3/5) <kF^2(r)> over the density: integrate numerically.
  double num = 0, den = 0;
  for (int i = 0; i < 4000; ++i) {
    const double r = (i + 0.5) * d.RMax() / 4000;
    const double w = r * r * d.Rho(r);
    const double kF = d.KFermi(pdg::kNeutron, r);
    num += w * 0.6 * kF * kF;
    den += w;
  }
  CHECK(sumK2 / N == Approx(num / den).epsilon(0.01));
  // Mean removal energy: S_n + 0.4 <T_F>, about 27-29 MeV for 12C.
  const double meanER = sumER / N;
  CHECK(meanER > 0.025);
  CHECK(meanER < 0.031);
  // Energy gap Q - Q_LFG = 16.827 MeV + (E_F^n - E_F^p) ~ 18.1 MeV, and V_out = T_F + S_n - gap.
  for (double r : {0.0, 1.0, 2.0, 3.0}) {
    const double gap = LocalFermiGas::EnergyGap(C, d, pdg::kNeutron, pdg::kProton, r);
    CHECK(gap == Approx(0.016827 + 0.00129).margin(5e-5));
    InitialState s;
    s.nHit = 1;
    s.hit[0].pdg = pdg::kNeutron;
    s.position = {0, 0, r};
    s.hasPosition = true;
    const double kF = d.KFermi(pdg::kProton, r);
    const double TF = std::sqrt(kF * kF + 0.938272 * 0.938272) - 0.938272;
    CHECK(lfg.OutgoingPotential(C, s, pdg::kProton) == Approx(TF + 0.0006).margin(1e-4));
    CHECK(lfg.PauliBlocked(C, pdg::kProton, {0, 0, 0.99 * kF}, s));
    CHECK_FALSE(lfg.PauliBlocked(C, pdg::kProton, {0, 0, 1.01 * kF}, s));
  }
  // Same-species transition (electron scattering, NC): no gap.
  CHECK(LocalFermiGas::EnergyGap(C, d, pdg::kProton, pdg::kProton, 1.0) == 0.0);
}

TEST_CASE("Lindhard function limits and imaginary part", "[lfg][rpa]") {
  const double M = 0.938919, kF = 0.25;
  const double rho = 2 * kF * kF * kF / (3 * constants::kPi * constants::kPi);
  // Static long-wavelength limit: -g M kF / (2 pi^2), g = 4.
  const auto u0 = ValenciaRPA::LindhardN(1e-7, 1e-3, kF, M);
  CHECK(u0.real() == Approx(-2 * M * kF / (constants::kPi * constants::kPi)).epsilon(1e-3));
  // High frequency: rho [1/(q0 - e) - 1/(q0 + e)], e = q^2/2M.
  const double q = 0.1, q0 = 0.5, e = q * q / (2 * M);
  CHECK(ValenciaRPA::LindhardN(q0, q, kF, M).real() == Approx(rho * (1 / (q0 - e) - 1 / (q0 + e))).epsilon(2e-3));
  // Im U_N(q0) = -pi rho < (1 - n(p+q)) delta(q0 - (E(p+q) - E(p))) >_{Fermi sphere}, by Monte Carlo.
  Random rng(5);
  const double qq = 0.3;
  const int nb = 40;
  const double wmax = 0.12;
  std::vector<double> hist(nb, 0.0);
  const int N = 2000000;
  for (int i = 0; i < N; ++i) {
    ThreeVector p;
    do p = RandomVec(rng, kF);
    while (p.Mag2() > kF * kF);
    const ThreeVector pq = p + ThreeVector{0, 0, qq};
    if (pq.Mag2() < kF * kF) continue;
    const double w = (pq.Mag2() - p.Mag2()) / (2 * M);
    const int b = static_cast<int>(w / wmax * nb);
    if (b >= 0 && b < nb) hist[b] += 1;
  }
  for (int b : {5, 10, 15, 20}) {
    const double w = (b + 0.5) * wmax / nb;
    const double mc = -constants::kPi * rho * hist[b] / N / (wmax / nb);
    CHECK(ValenciaRPA::LindhardN(w, qq, kF, M).imag() == Approx(mc).epsilon(0.03));
  }
  // Complex-energy continuation: tends to the real-axis function (inside and outside the support).
  for (double w : {0.02, 0.06, 0.1, 0.2, 0.4}) {
    const auto u = ValenciaRPA::LindhardN(w, qq, kF, M), uc = ValenciaRPA::LindhardN(w, qq, kF, M, 1e-6);
    CHECK(uc.real() == Approx(u.real()).epsilon(1e-3).margin(1e-6));
    CHECK(uc.imag() == Approx(u.imag()).epsilon(1e-3).margin(1e-6));
  }
}

TEST_CASE("RPA factors: limits", "[rpa]") {
  ValenciaRPA rpa;
  // No medium: no correction.
  const auto f0 = rpa.Evaluate(0.05, 0.3, 0.0, 0.0);
  CHECK(f0.CN == 1.0);
  CHECK(f0.CL == 1.0);
  CHECK(f0.CT == 1.0);
  // Low q at nuclear-matter density: strong quenching; large q: corrections fade away.
  const auto lo = rpa.Evaluate(0.02, 0.1, 0.16, 1.0);
  CHECK(lo.CT < 0.6);
  CHECK(lo.CN < 0.8);
  // Beyond the validity range of the non-relativistic interaction (|q| >= q_zero = 1 GeV): no correction;
  // between q_full and q_zero the corrections fade linearly.
  const auto hi = rpa.Evaluate(0.8, 1.5, 0.16, 1.0);
  CHECK(hi.CT == 1.0);
  CHECK(hi.CL == 1.0);
  CHECK(hi.CN == 1.0);
  ValenciaRPA::Parameters p;
  p.qZero = 10;
  p.qFull = 10;
  const ValenciaRPA full(p);
  const auto mid = rpa.Evaluate(0.34, 0.8, 0.16, 1.0), midFull = full.Evaluate(0.34, 0.8, 0.16, 1.0);
  CHECK(mid.CT - 1 == Approx(0.5 * (midFull.CT - 1)).epsilon(1e-9));
}

TEST_CASE("RPA block contraction matches brute-force tensor decomposition", "[rpa]") {
  static const test::Gammas G;
  const std::complex<double> I(0, 1);
  Random rng(17);
  const double M = constants::kMassNucleon, m = constants::kMassMuon;
  // Hadron tensor H^{mu nu} by explicit traces for given form factors (axialSign = -1 as in test_qe_current).
  auto hadron = [&](const FourVector& p, const FourVector& pp, double F1, double F2, double FA, double FP,
                    std::complex<double> H[4][4]) {
    const FourVector q = pp - p;
    test::Mat Gam[4], GamBar[4];
    for (int mu = 0; mu < 4; ++mu) {
      test::Mat sigq;
      for (int al = 0; al < 4; ++al) {
        const test::Mat sig = (G.g[mu] * G.g[al] - G.g[al] * G.g[mu]) * (0.5 * I);
        sigq = sigq + sig * (test::Metric(al) * test::Comp(q, al));
      }
      Gam[mu] = G.g[mu] * F1 + sigq * (I * F2 / (2 * M)) + G.g[mu] * G.g5 * (-FA) + G.g5 * (-test::Comp(q, mu) * FP / M);
      GamBar[mu] = G.Bar(Gam[mu]);
    }
    const test::Mat Pp = G.Slash(pp) + G.one * M, Pi = G.Slash(p) + G.one * M;
    for (int mu = 0; mu < 4; ++mu)
      for (int nu = 0; nu < 4; ++nu) H[mu][nu] = (Pp * Gam[mu] * Pi * GamBar[nu]).Trace();
  };
  auto block = [](int a, int b) {  // 0:00 1:0L 2:LL 3:TT 4:0T 5:LT (q along z)
    const bool t0 = a == 0, t1 = b == 0, l0 = a == 3, l1 = b == 3;
    if (t0 && t1) return 0;
    if ((t0 && l1) || (l0 && t1)) return 1;
    if (l0 && l1) return 2;
    if (t0 || t1) return 4;
    if (l0 || l1) return 5;
    return 3;
  };
  for (int trial = 0; trial < 20; ++trial) {
    // Kinematics with q along z: k, k' random, p random, p' = p + q on shell for a q0 that follows.
    const FourVector k = OnShell(0, ThreeVector{0.3 * (2 * rng.Uniform() - 1), 0.3 * (2 * rng.Uniform() - 1), 2.0});
    const double qz = 0.2 + 0.8 * rng.Uniform();
    const ThreeVector kpv = k.Vect() - ThreeVector{0, 0, qz};
    const FourVector kp = OnShell(m, kpv);
    const FourVector p = OnShell(M, RandomVec(rng, 0.2));
    const FourVector pp = OnShell(M, p.Vect() + ThreeVector{0, 0, qz});
    const WeakFormFactors ff{rng.Uniform(0.2, 1.0), rng.Uniform(0.5, 3.0), rng.Uniform(0.3, 1.2), rng.Uniform(1, 20)};
    const ValenciaRPA::Factors fac{0.7, 0.45, 0.3};
    for (int h : {+1, -1}) {
      // Lepton tensor by explicit traces.
      const test::Mat PL = G.one - G.g5;
      const test::Mat Kin = h > 0 ? G.Slash(k) : G.Slash(kp) - G.one * m;
      const test::Mat Kout = h > 0 ? G.Slash(kp) + G.one * m : G.Slash(k);
      std::complex<double> L[4][4];
      for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu) L[mu][nu] = (Kout * G.g[mu] * PL * Kin * G.g[nu] * PL).Trace();
      // Structures by polarization of the bilinear form-factor dependence.
      std::complex<double> H11[4][4], H22[4][4], HAA[4][4], HPP[4][4], H12s[4][4], HAPs[4][4], H1As[4][4],
          H2As[4][4], H1Ps[4][4], H2Ps[4][4];
      hadron(p, pp, ff.F1, 0, 0, 0, H11);
      hadron(p, pp, 0, ff.F2, 0, 0, H22);
      hadron(p, pp, 0, 0, ff.FA, 0, HAA);
      hadron(p, pp, 0, 0, 0, ff.FP, HPP);
      hadron(p, pp, ff.F1, ff.F2, 0, 0, H12s);
      hadron(p, pp, 0, 0, ff.FA, ff.FP, HAPs);
      hadron(p, pp, ff.F1, 0, ff.FA, 0, H1As);
      hadron(p, pp, 0, ff.F2, ff.FA, 0, H2As);
      hadron(p, pp, ff.F1, 0, 0, ff.FP, H1Ps);
      hadron(p, pp, 0, ff.F2, 0, ff.FP, H2Ps);
      const double CN = fac.CN, CT = fac.CT, CL = fac.CL;
      const double f11[6] = {CN, CN, CN, 1, 1, 1}, f12[6] = {CN, CN, CN, CT, 1, 1}, f22[6] = {1, 1, 1, CT, 1, 1};
      const double fAA[6] = {1, CL, CL, CT, 1, 1}, fP[6] = {CL, CL, CL, 1, 1, 1}, fVA[6] = {1, 1, 1, CT, 1, 1};
      std::complex<double> tot = 0, totNoRPA = 0;
      for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu) {
          const int b = block(mu, nu);
          const std::complex<double> h12 = H12s[mu][nu] - H11[mu][nu] - H22[mu][nu];
          const std::complex<double> hAP = HAPs[mu][nu] - HAA[mu][nu] - HPP[mu][nu];
          const std::complex<double> hVA = (H1As[mu][nu] - H11[mu][nu] - HAA[mu][nu]) + (H2As[mu][nu] - H22[mu][nu] - HAA[mu][nu]);
          const std::complex<double> hVP = (H1Ps[mu][nu] - H11[mu][nu] - HPP[mu][nu]) + (H2Ps[mu][nu] - H22[mu][nu] - HPP[mu][nu]);
          const std::complex<double> hr = f11[b] * H11[mu][nu] + f12[b] * h12 + f22[b] * H22[mu][nu] + fAA[b] * HAA[mu][nu] +
                                          fP[b] * (hAP + HPP[mu][nu]) + fVA[b] * hVA + hVP;
          const std::complex<double> hall = H11[mu][nu] + h12 + H22[mu][nu] + HAA[mu][nu] + hAP + HPP[mu][nu] + hVA + hVP;
          const double gg = test::Metric(mu) * test::Metric(nu);
          tot += gg * L[mu][nu] * hr;
          totNoRPA += gg * L[mu][nu] * hall;
        }
      // Random rotation of the whole event: the contraction must not depend on the orientation.
      const double a = 2 * constants::kPi * rng.Uniform(), c = std::cos(a), s = std::sin(a);
      auto rot = [&](const FourVector& v) {  // rotation about the x axis
        return FourVector(v.e, v.px, c * v.py - s * v.pz, s * v.py + c * v.pz);
      };
      const double code = ValenciaRPA::Contract(rot(k), rot(kp), rot(p), rot(pp), ff, M, h, fac);
      const double codeNoRPA = ValenciaRPA::Contract(rot(k), rot(kp), rot(p), rot(pp), ff, M, h, ValenciaRPA::Factors{});
      CHECK(codeNoRPA == Approx(ContractQE(k, kp, p, pp, ff, M, h)).epsilon(1e-9));
      CHECK(codeNoRPA == Approx(totNoRPA.real()).epsilon(1e-9));
      CHECK(code == Approx(tot.real()).epsilon(1e-8));
    }
  }
}

TEST_CASE("CCQE events with LFG / SF + RPA and the cascade conserve four-momentum", "[lfg][generator]") {
  for (const char* nm : {"LFG", "SF"}) {
  YAML::Node c = YAML::Load(R"(
run: {seed: 7, batch_size: 1000, threads: 4}
beam: {pdg: 14, flux: {type: mono, energy: 3.0}}
target: "1000060120:1"
channels:
  - type: CCQE
    nuclear_model: {type: XX}
    form_factors: {vector: BBBA05, axial: dipole, MA: 1.03}
    rpa: true
fsi: {type: cascade}
xsec_table: {nodes: 10, samples: 20000, cache: none}
)");
  c["channels"][0]["nuclear_model"]["type"] = nm;
  struct W : OutputWriter {
    std::vector<Event> ev;
    void Begin(const RunInfo&) override {}
    void Write(const Event& e) override { ev.push_back(e); }
    void End(const RunSummary&) override {}
  } w;
  Generator gen(c);
  gen.Run(3000, {&w});
  REQUIRE(w.ev.size() == 3000);
  for (const auto& e : w.ev) {
    FourVector in, outPre, outPost;
    ThreeVector vtx;
    for (const auto& p : e.particles) {
      if (p.status == kStatusBeam || p.status == kStatusTarget) in = in + p.p;
      if (p.status == kStatusPreFSI || p.status == kStatusRemnant || (p.status == kStatusFinal && std::abs(p.pdg) == 13))
        outPre = outPre + p.p;
      if (p.status == kStatusFinal || p.status == kStatusRemnantPostFSI) outPost = outPost + p.p;
      if (p.status == kStatusPreFSI) vtx = p.x;
    }
    REQUIRE(std::abs(in.e - outPre.e) < 1e-9);
    REQUIRE((in.Vect() - outPre.Vect()).Mag() < 1e-9);
    REQUIRE(std::abs(in.e - outPost.e) < 1e-9);
    REQUIRE((in.Vect() - outPost.Vect()).Mag() < 1e-9);
    // The cascade starts from the stage-2 interaction point.
    REQUIRE(e.initial.hasPosition);
    REQUIRE((vtx - e.initial.position).Mag() < 1e-12);
  }
  }
}

TEST_CASE("12C spectral function: shells, SRC pairs and energy balance", "[sf]") {
  const Nucleus& C = NuclearDatabase::Instance().Get(1000060120);
  ShellModelSpectralFunction sf;
  const auto& par = sf.Params();
  Random rng(17);
  HitNucleonRequest req;
  req.species[0] = pdg::kNeutron;
  const int N = 300000;
  int nSRC = 0, nNP = 0;
  double sumERmf = 0, sumK2mf = 0;
  int nMF = 0;
  for (int i = 0; i < N; ++i) {
    InitialState s;
    sf.Sample(C, req, s, rng);
    REQUIRE(s.hasPosition);
    FourVector tot = s.hit[0].Momentum() + s.remnant;
    if (s.hasSpectator) {
      ++nSRC;
      tot = tot + s.spectator.Momentum();
      REQUIRE(s.hit[0].p.Mag() >= par.kSRC - 1e-12);
      REQUIRE(s.hit[0].p.Mag() <= par.kMax + 1e-12);
      if (s.spectator.pdg == pdg::kProton) ++nNP;
      // A-2 residual in its ground state: 10B for an np pair.
      if (s.spectator.pdg == pdg::kProton) REQUIRE(s.remnantPdg == pdg::IonCode(5, 10));
    } else {
      ++nMF;
      sumERmf += s.hit[0].removalEnergy;
      sumK2mf += s.hit[0].p.Mag2();
      REQUIRE(s.hit[0].removalEnergy >= C.separationNeutron - 1e-9);
      REQUIRE(s.remnantPdg == pdg::IonCode(6, 11));
    }
    // Nucleus at rest: the hit nucleon, the partner and the residual carry M_A and zero momentum.
    REQUIRE(std::abs(tot.e - C.mass) < 1e-9);
    REQUIRE(tot.Vect().Mag() < 1e-9);
  }
  CHECK(static_cast<double>(nSRC) / N == Approx(1 - par.mfFraction).margin(0.005));
  CHECK(static_cast<double>(nNP) / nSRC == Approx(par.npFraction).margin(0.005));
  // Mean field: <k^2> of the HO 1s (3/2) and 1p (5/2) orbitals, weights 2/6 and 4/6, in units of (hbar c / b)^2.
  const double b = NuclearDensity::For(C).B();
  const double k2HO = (2.0 / 6 * 1.5 + 4.0 / 6 * 2.5) * (kHbarC / b) * (kHbarC / b);
  CHECK(sumK2mf / nMF == Approx(k2HO).epsilon(0.01));
  // Mean neutron removal energy of the mean field: 1/3 1s + 2/3 1p centroids shifted by S_n - S_p, plus
  // the upward skew of the Lorentzians whose width grows with E (tails cut at E_max): 0-6 MeV above.
  const double dS = C.separationNeutron - C.separationProton;
  const double naive = (par.e1s + 2 * par.e1p) / 3 + dS;
  CHECK(sumERmf / nMF > naive);
  CHECK(sumERmf / nMF < naive + 0.006);
  // Medians of the shell distributions sit at the centroids.
  CHECK(sf.SampleRemovalEnergy(1, 0.5) == Approx(par.e1p).margin(0.001));
  CHECK(sf.SampleRemovalEnergy(0, 0.5) == Approx(par.e1s).margin(0.004));
}
