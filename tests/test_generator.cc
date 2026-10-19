// End-to-end tests of the orchestration layer.
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#include "core/PDG.h"
#include "orchestration/Generator.h"

using namespace nuflow;
using Catch::Approx;

namespace {

YAML::Node BaseConfig(const std::string& flux, const std::string& target, int nu = 14, int threads = 4) {
  YAML::Node c = YAML::Load(R"(
run: {seed: 99, batch_size: 1000}
channels:
  - type: CCQE
    nuclear_model: {type: RFG}
    form_factors: {vector: BBBA05, axial: dipole, MA: 1.03}
fsi: {type: none}
xsec_table: {nodes: 40, samples: 40000, cache: none}
)");
  c["run"]["threads"] = threads;
  c["beam"]["pdg"] = nu;
  c["beam"]["flux"] = YAML::Load(flux);
  c["target"] = target;
  return c;
}

// Collects events in memory.
struct MemoryWriter : OutputWriter {
  std::vector<Event> events;
  void Begin(const RunInfo&) override {}
  void Write(const Event& e) override { events.push_back(e); }
  void End(const RunSummary&) override {}
};

double Chi2PerDof(const std::vector<double>& obs, const std::vector<double>& exp, int* ndf) {
  double chi2 = 0;
  *ndf = 0;
  for (size_t i = 0; i < obs.size(); ++i) {
    if (exp[i] < 10) continue;
    chi2 += (obs[i] - exp[i]) * (obs[i] - exp[i]) / exp[i];
    ++*ndf;
  }
  return chi2 / *ndf;
}

}  // namespace

TEST_CASE("Events conserve four-momentum and have the right particles", "[generator]") {
  for (int nu : {14, -14}) {
    Generator gen(BaseConfig("{type: mono, energy: 1.5}", "1000060120:1,1000010010:2", nu));
    MemoryWriter w;
    gen.Run(5000, {&w});
    REQUIRE(w.events.size() == 5000);
    int nHydrogen = 0;
    for (const auto& ev : w.events) {
      FourVector in, out;
      int nLepton = 0;
      for (const auto& p : ev.particles) {
        if (p.status == kStatusBeam || p.status == kStatusTarget) in = in + p.p;
        if (p.status == kStatusFinal || p.status == kStatusRemnant) out = out + p.p;
        nLepton += p.status == kStatusFinal && std::abs(p.pdg) == 13;
      }
      REQUIRE(nLepton == 1);
      REQUIRE(std::abs(in.e - out.e) < 1e-9);
      REQUIRE((in.Vect() - out.Vect()).Mag() < 1e-9);
      nHydrogen += ev.targetPdg == pdg::kProton;
      REQUIRE(ev.enu == Approx(1.5));
    }
    if (nu > 0) CHECK(nHydrogen == 0);  // no QE on hydrogen for neutrinos
    else CHECK(nHydrogen > 0);
  }
}

TEST_CASE("Output is independent of the number of threads", "[generator]") {
  MemoryWriter a, b;
  Generator(BaseConfig("{type: mono, energy: 2.0}", "1000060120:1", 14, 1)).Run(3000, {&a});
  Generator(BaseConfig("{type: mono, energy: 2.0}", "1000060120:1", 14, 7)).Run(3000, {&b});
  REQUIRE(a.events.size() == b.events.size());
  for (size_t i = 0; i < a.events.size(); ++i) {
    REQUIRE(a.events[i].Q2 == b.events[i].Q2);
    REQUIRE(a.events[i].particles.size() == b.events[i].particles.size());
  }
}

TEST_CASE("Unweighted Q2 distribution matches the weighted trial distribution", "[generator]") {
  const double E = 3.0;
  Generator gen(BaseConfig("{type: mono, energy: 3.0}", "1000060120:1"));
  const int nb = 30;
  const double q2max = 1.5;
  std::vector<double> unw(nb, 0), wtd(nb, 0);
  MemoryWriter w;
  gen.Run(100000, {&w});
  for (const auto& ev : w.events)
    if (ev.Q2 < q2max) unw[static_cast<int>(ev.Q2 / q2max * nb)] += 1;

  const Process& proc = gen.Selector().Get(0);
  Random rng(1234);
  Event ev;
  double wsum = 0;
  for (int i = 0; i < 400000; ++i) {
    const double wt = proc.channel->Trial(14, E, proc.nucleus, ev, rng);
    if (wt <= 0) continue;
    wsum += wt;
    if (ev.Q2 < q2max) wtd[static_cast<int>(ev.Q2 / q2max * nb)] += wt;
  }
  for (auto& x : wtd) x *= w.events.size() / wsum;
  int ndf;
  const double chi2 = Chi2PerDof(unw, wtd, &ndf);
  INFO("chi2/ndf = " << chi2 << " ndf = " << ndf);
  CHECK(chi2 < 2.0);
}

TEST_CASE("Event energy spectrum follows flux x cross section", "[generator]") {
  const std::string flux =
      "{type: histogram, format: root, file: data/flux/minerva_le_2016/minerva_flux.root, histogram: numu_fhc,"
      " density: per_bin, scale: 1.0e-6, emax: 20}";
  Generator gen(BaseConfig(flux, "1000060120:1,1000010010:1"));
  MemoryWriter w;
  gen.Run(200000, {&w});
  const auto f = MakeFlux(YAML::Load(flux));
  const auto& h = dynamic_cast<const HistogramFlux&>(*f);
  std::vector<double> obs(h.NBins(), 0), exp(h.NBins(), 0);
  for (const auto& ev : w.events) {
    const auto it = std::upper_bound(h.Edges().begin(), h.Edges().end(), ev.enu);
    obs[std::min<size_t>(it - h.Edges().begin() - 1, h.NBins() - 1)] += 1;
  }
  double tot = 0;
  for (size_t i = 0; i < h.NBins(); ++i) {
    // fine numerical integral of sigma over the bin
    const double lo = h.Edges()[i], hi = h.Edges()[i + 1];
    double s = 0;
    for (int j = 0; j < 200; ++j) s += gen.Selector().TotalSigma(lo + (hi - lo) * (j + 0.5) / 200) * (hi - lo) / 200;
    exp[i] = h.Density()[i] * s;
    tot += exp[i];
  }
  for (auto& x : exp) x *= w.events.size() / tot;
  int ndf;
  const double chi2 = Chi2PerDof(obs, exp, &ndf);
  INFO("chi2/ndf = " << chi2 << " ndf = " << ndf);
  CHECK(chi2 < 2.0);
  // Flux-averaged cross section is consistent with the same integral.
  double fluxInt = 0, num = 0;
  for (size_t i = 0; i < h.NBins(); ++i) fluxInt += h.Density()[i] * (h.Edges()[i + 1] - h.Edges()[i]);
  for (size_t i = 0; i < h.NBins(); ++i) num += exp[i] * tot / w.events.size();
  CHECK(gen.Sampler().FluxAveragedSigma() == Approx(num / fluxInt).epsilon(1e-3));
}
