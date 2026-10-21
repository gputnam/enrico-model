// Validation of the 2p2h meson-exchange-current model (theory/mec.md):
//  1. label symmetries of the isospin-stripped currents used to derive the charge-channel formulas
//     (RS17 Eqs. 57-63): J1(2'1';21) = J2(1'2';12), J3(2'1';21) = -J3(1'2';12);
//  2. the vector Delta MEC (C3V vertex) is conserved, Q_mu J^mu = 0;
//  3. the response functions of 12C reproduce the published RS17 values (Figs. 7, 8, 11) and the
//     pp/np ratio of Ruiz Simo et al. PLB 762 (2016) 124;
//  4. events: four-momentum and charge conservation, two nucleons above the Fermi surface, process id 300;
//  5. sigma(E) saturates and sigma(nu) > sigma(nubar) (constructive V-A interference for neutrinos).
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "core/Constants.h"
#include "core/Dirac.h"
#include "core/PDG.h"
#include "core/Random.h"
#include "kinematics/MECCurrent.h"
#include "orchestration/Generator.h"
#include "orchestration/OutputWriter.h"

using namespace nuflow;
using Catch::Approx;
using cplx = std::complex<double>;

namespace {

NucleonFormFactors HNVNucleonFF() {
  return NucleonFormFactors::FromConfig(YAML::Load("{vector: dipole, axial: dipole, gA: 1.26, MA: 1.05}"));
}

FourVector OnShell(const ThreeVector& p) {
  const double M = constants::kMassNucleon;
  return {std::sqrt(p.Mag2() + M * M), p};
}

ThreeVector RandomVec(Random& r, double scale) {
  return {scale * (2 * r.Uniform() - 1), scale * (2 * r.Uniform() - 1), scale * (2 * r.Uniform() - 1)};
}

int SwapSpins(int s) {
  const int sf1 = s & 1, sf2 = (s >> 1) & 1, si1 = (s >> 2) & 1, si2 = (s >> 3) & 1;
  return sf2 | (sf1 << 1) | (si2 << 2) | (si1 << 3);
}

YAML::Node MecConfig(const std::string& target, double E, int nu) {
  YAML::Node c = YAML::Load(R"(
run: {seed: 11, batch_size: 500, threads: 4}
channels:
  - type: CCMEC
    nuclear_model: {type: RFG}
fsi: {type: none}
xsec_table: {samples: 100000, cache: none}
)");
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

struct Resp {
  double CC, T, Tp;
};
Resp Responses(const MECCurrent& cur, MECChannel ch, double q, double w, long n, bool direct = false) {
  Random rng(17, static_cast<uint64_t>(1000 * w + 7 * q));
  cplx W[4][4];
  MECHadronTensor(cur, ch, q, w - 0.040, 0.228, 6, n, rng, W, direct);  // RS17 Fig. 7: kF = 228 MeV, eps = 40 MeV
  return {W[0][0].real(), (W[1][1] + W[2][2]).real(), (cplx(0, -0.5) * (W[1][2] - W[2][1])).real()};
}

}  // namespace

TEST_CASE("MEC currents: label symmetries behind the charge-channel formulas", "[mec]") {
  const MECCurrent cur(MECParameters{}, HNVNucleonFF());
  Random rng(3);
  for (int trial = 0; trial < 5; ++trial) {
    const FourVector h1 = OnShell(RandomVec(rng, 0.13)), h2 = OnShell(RandomVec(rng, 0.13));
    const FourVector p1 = OnShell(RandomVec(rng, 0.6) + ThreeVector{0, 0, 0.3});
    const FourVector p2 = OnShell(RandomVec(rng, 0.6));
    const FourVector Q = p1 + p2 - h1 - h2;
    MECCurrent::Currents a, b;
    cur.ComputeCurrents(Q, p1, h1, p2, h2, a);
    cur.ComputeCurrents(Q, p2, h2, p1, h1, b);  // full 1 <-> 2 relabelling
    double scale = 0;
    for (int s = 0; s < 16; ++s)
      for (int mu = 0; mu < 4; ++mu) scale = std::max(scale, std::abs(a.J3[s][mu]) + std::abs(a.J2[s][mu]));
    for (int s = 0; s < 16; ++s)
      for (int mu = 0; mu < 4; ++mu) {
        CHECK(std::abs(b.J1[SwapSpins(s)][mu] - a.J2[s][mu]) < 1e-10 * scale);
        CHECK(std::abs(b.J3[SwapSpins(s)][mu] + a.J3[s][mu]) < 1e-10 * scale);
      }
  }
}

TEST_CASE("Vector Delta MEC is conserved", "[mec]") {
  MECParameters par;
  par.axialScale = 0;
  par.backgroundScale = 0;  // Delta terms only (C3V vertex is transverse to Q)
  const MECCurrent cur(par, HNVNucleonFF());
  Random rng(5);
  for (int trial = 0; trial < 5; ++trial) {
    const FourVector h1 = OnShell(RandomVec(rng, 0.13)), h2 = OnShell(RandomVec(rng, 0.13));
    const FourVector p1 = OnShell(RandomVec(rng, 0.6) + ThreeVector{0, 0, 0.4}), p2 = OnShell(RandomVec(rng, 0.5));
    const FourVector Q = p1 + p2 - h1 - h2;
    MECCurrent::Currents c;
    cur.ComputeCurrents(Q, p1, h1, p2, h2, c);
    for (int s = 0; s < 16; ++s) {
      cplx d1 = 0, d3 = 0;
      double n1 = 0, n3 = 0;
      for (int mu = 0; mu < 4; ++mu) {
        const double qmu = dirac::Metric(mu) * dirac::Comp(Q, mu);
        d1 += qmu * c.J1[s][mu];
        d3 += qmu * c.J3[s][mu];
        n1 += std::abs(c.J1[s][mu]) * std::abs(dirac::Comp(Q, mu));
        n3 += std::abs(c.J3[s][mu]) * std::abs(dirac::Comp(Q, mu));
      }
      CHECK(std::abs(d1) <= 1e-10 * (n1 + 1e-30));
      CHECK(std::abs(d3) <= 1e-10 * (n3 + 1e-30));
    }
  }
}

TEST_CASE("12C 2p2h responses vs Ruiz Simo et al. (2017)", "[mec]") {
  const MECCurrent cur(MECParameters{}, HNVNucleonFF());
  const long n = 6000;
  // RS17 Fig. 8 (q = 800 MeV/c, omega = 500 MeV): R^T ~ 15, R^T' ~ 7.5, R^CC ~ 2.6 GeV^-1.
  const Resp pp = Responses(cur, MECChannel::kNuNPtoPP, 0.8, 0.5, n);
  const Resp np = Responses(cur, MECChannel::kNuNNtoNP, 0.8, 0.5, n);
  const double T = pp.T + np.T, Tp = pp.Tp + np.Tp, CC = pp.CC + np.CC;
  INFO("q=800, w=500: R_T " << T << " R_T' " << Tp << " R_CC " << CC << " pp/np(T) " << pp.T / np.T);
  CHECK(T == Approx(15.0).epsilon(0.2));
  CHECK(Tp == Approx(7.5).epsilon(0.2));
  CHECK(CC == Approx(2.6).epsilon(0.3));
  // PLB 762 (2016) 124: pp emission dominates the transverse response by "a factor 6 or less".
  CHECK(pp.T / np.T > 3.5);
  CHECK(pp.T / np.T < 8.0);
  // RS17 Fig. 11 (q = 400 MeV/c, omega ~ 300 MeV): direct only ~71, direct + exchange ~51 GeV^-1.
  const Resp full = Responses(cur, MECChannel::kNuNPtoPP, 0.4, 0.3, n);
  const Resp fullnp = Responses(cur, MECChannel::kNuNNtoNP, 0.4, 0.3, n);
  const Resp dir = Responses(cur, MECChannel::kNuNPtoPP, 0.4, 0.3, n, true);
  const Resp dirnp = Responses(cur, MECChannel::kNuNNtoNP, 0.4, 0.3, n, true);
  const double ratio = (dir.T + dirnp.T) / (full.T + fullnp.T);
  INFO("q=400, w=300: R_T " << full.T + fullnp.T << ", direct/total " << ratio);
  CHECK(full.T + fullnp.T == Approx(51).epsilon(0.2));
  CHECK(ratio == Approx(71.0 / 51.0).epsilon(0.1));
}

TEST_CASE("CCMEC events: conservation laws and two-nucleon final states", "[mec]") {
  for (int nu : {14, -14}) {
    Generator gen(MecConfig("1000060120:1", 3.0, nu));
    MemWriter w;
    gen.Run(2000, {&w});
    const double kF = NuclearDatabase::Instance().Get(1000060120).kFermiProton;
    int nPP = 0;
    for (const auto& ev : w.events) {
      FourVector in, out;
      int nLep = 0, nP = 0, nN = 0, nOther = 0, nHit = 0;
      for (const auto& p : ev.particles) {
        if (p.status == kStatusBeam || p.status == kStatusTarget) in = in + p.p;
        if (p.status == kStatusFinal || p.status == kStatusRemnant) out = out + p.p;
        if (p.status == kStatusStruckNucleon) ++nHit;
        if (p.status != kStatusFinal) continue;
        if (std::abs(p.pdg) == 13) {
          ++nLep;
          CHECK(p.pdg == (nu > 0 ? 13 : -13));
        } else if (p.pdg == 2212 || p.pdg == 2112) {
          (p.pdg == 2212 ? nP : nN)++;
          CHECK(p.p.P() > kF - 1e-9);
        } else {
          ++nOther;
        }
      }
      CHECK(nLep == 1);
      CHECK(nHit == 2);
      CHECK(nP + nN == 2);
      CHECK(nOther == 0);
      if (nu > 0) CHECK(nP >= 1);   // pp or np
      else CHECK(nN >= 1);          // np or nn
      if (nP == 2) ++nPP;
      CHECK(std::abs(in.e - out.e) < 1e-7);
      CHECK((in.Vect() - out.Vect()).Mag() < 1e-7);
      CHECK(ev.channel == 300);
    }
    if (nu > 0) {
      INFO("pp fraction " << double(nPP) / w.events.size());
      CHECK(double(nPP) / w.events.size() > 0.6);  // np-pair (pp final state) dominance
    }
  }
}

TEST_CASE("CCMEC cross section: saturation with energy and nu > nubar", "[mec]") {
  Generator g3(MecConfig("1000060120:1", 3.0, 14));
  Generator g10(MecConfig("1000060120:1", 10.0, 14));
  Generator gb(MecConfig("1000060120:1", 3.0, -14));
  const double s3 = g3.Selector().Get(0).table.sigma[0], s10 = g10.Selector().Get(0).table.sigma[0];
  const double sb = gb.Selector().Get(0).table.sigma[0];
  INFO("sigma_2p2h(C12) [1e-38 cm2]: nu 3 GeV " << s3 << ", nu 10 GeV " << s10 << ", nubar 3 GeV " << sb);
  CHECK(s3 > 1.0);
  CHECK(s3 < 2.5);
  CHECK(s10 / s3 == Approx(1.0).epsilon(0.15));
  CHECK(sb < s3);
}
