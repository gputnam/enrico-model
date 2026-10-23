// nuflow-qe-table: reproduces the CCQE validation table of theory/quasielastic.md.
// For each energy prints, per struck nucleon [1e-38 cm^2]:
//   the Llewellyn Smith dsigma/dQ2 integrated over Q2 (free nucleon at rest, average nucleon mass),
//   the nuflow Monte Carlo for a free nucleon, the nuflow Monte Carlo for 12C (RFG), and C/free,
// for nu_mu n and nubar_mu p, with the dipole (M_A = 1.03 GeV), the deuterium z-expansion
// (Meyer et al. 2016) and the LQCD-average z-expansion (Meyer, arXiv:2601.02676) axial form factors.
//
//   nuflow-qe-table [samples per point, default 1e6]

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "core/Constants.h"
#include "core/Nucleus.h"
#include "core/PDG.h"
#include "initial_state/FreeNucleon.h"
#include "initial_state/RelativisticFermiGas.h"
#include "kinematics/QECurrent.h"
#include "kinematics/QuasiElasticKinematics.h"

using namespace nuflow;

static double IntegratedLS(double E, const NucleonFormFactors& ff, int h) {
  using namespace constants;
  const double M = kMassNucleon, m = kMassMuon, s = M * M + 2 * M * E, rs = std::sqrt(s);
  const double kc = (s - M * M) / (2 * rs), El = (s + m * m - M * M) / (2 * rs);
  if (El <= m) return 0;
  const double pl = std::sqrt(El * El - m * m), a = -m * m + 2 * kc * (El - pl), b = -m * m + 2 * kc * (El + pl);
  const int n = 4000;
  double sum = 0;
  for (int i = 0; i <= n; ++i) {
    const double Q2 = a + (b - a) * i / n, w = (i == 0 || i == n) ? 1 : (i % 2 ? 4 : 2);
    sum += w * LlewellynSmithDSigmaDQ2(E, Q2, m, M, ff.Evaluate(Q2), h);
  }
  return sum * (b - a) / (3 * n) * kGeV2To1e38Cm2;
}

static double MonteCarlo(int nu, double E, const Nucleus& A, const NuclearModel& model,
                         const QuasiElasticKinematics& kin, long N) {
  Random rng(11, static_cast<uint64_t>(E * 1e4));
  Event ev;
  HitNucleonRequest req;
  req.species[0] = nu > 0 ? pdg::kNeutron : pdg::kProton;
  double sum = 0;
  for (long i = 0; i < N; ++i) {
    InitialState init;
    model.Sample(A, req, init, rng);
    sum += kin.Generate(nu, E, A, init, model, ev, rng);
  }
  return sum / N;
}

int main(int argc, char** argv) {
  const long N = argc > 1 ? std::atol(argv[1]) : 1000000;
  const auto& db = NuclearDatabase::Instance();
  FreeNucleon free;
  RelativisticFermiGas rfg;
  for (const std::string set : {"dipole", "zexp", "zexp_lqcd"}) {
    YAML::Node c;
    c["vector"] = "BBBA05";
    c["axial"] = set == "dipole" ? "dipole" : "zexp";
    c["MA"] = 1.03;
    if (set == "zexp_lqcd") {
      c["gA"] = 1.2754;
      c["tcut"] = 9 * 0.134 * 0.134;
      c["t0"] = -0.50;
      c["zexp_a"] = std::vector<double>{-1.70104640, 0.26324902};
    }
    const char* axial = set.c_str();
    const auto ff = NucleonFormFactors::FromConfig(c);
    QuasiElasticKinematics kin(ff);
    std::printf("# axial=%s\n# E   LS(nu n)  MC(nu n)  C12/neutron  ratio   LS(nubar p) MC(nubar p) C12/proton ratio\n",
                axial);
    for (double E : {0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0}) {
      const double l1 = IntegratedLS(E, ff, 1), m1 = MonteCarlo(14, E, db.Get(2112), free, kin, N);
      const double c1 = MonteCarlo(14, E, db.Get(1000060120), rfg, kin, N);
      const double l2 = IntegratedLS(E, ff, -1), m2 = MonteCarlo(-14, E, db.Get(2212), free, kin, N);
      const double c2 = MonteCarlo(-14, E, db.Get(1000060120), rfg, kin, N);
      std::printf("%5.2f  %.4f  %.4f  %.4f  %.3f   %.4f  %.4f  %.4f  %.3f\n", E, l1, m1, c1, c1 / m1, l2, m2, c2,
                  c2 / m2);
    }
  }
  return 0;
}
