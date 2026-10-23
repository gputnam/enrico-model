// nuflow-lfg-validate: validation tables of the local Fermi gas and the Valencia RPA
// (theory/local_fermi_gas.md, Sec. 6).
//
//   nuflow-lfg-validate factors              RPA factors C_N, C_T, C_L along the QE peak at rho_0
//   nuflow-lfg-validate momentum             12C struck-nucleon momentum distribution n(k), RFG vs LFG
//   nuflow-lfg-validate q2 E [N] [nu]        CCQE dsigma/dQ2 on 12C at E: RFG, LFG, LFG+RPA and ratios
//   nuflow-lfg-validate xsec [N] [nu]        CCQE sigma(E) per neutron (proton) on 12C: RFG, LFG, LFG+RPA
//   nuflow-lfg-validate eqe E theta [N]      12C(e,e') QE d2sigma/dOmega/domega [nb/sr/GeV] vs omega,
//                                            RFG and LFG (vector current, no RPA)

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "core/Constants.h"
#include "core/Nucleus.h"
#include "core/PDG.h"
#include "core/Random.h"
#include "fsi/NuclearDensity.h"
#include "initial_state/LocalFermiGas.h"
#include "initial_state/RelativisticFermiGas.h"
#include "initial_state/ShellModelSpectralFunction.h"
#include "kinematics/QECurrent.h"
#include "kinematics/QuasiElasticKinematics.h"
#include "kinematics/ValenciaRPA.h"

using namespace nuflow;
using constants::kPi;

namespace {

const Nucleus& Carbon() { return NuclearDatabase::Instance().Get(1000060120); }

NucleonFormFactors FF() {
  return NucleonFormFactors::FromConfig(YAML::Load("{vector: BBBA05, axial: dipole, MA: 1.03}"));
}

void Factors() {
  ValenciaRPA rpa;
  const double M = 0.938919, rho0 = 0.16;
  std::printf("# RPA factors at rho = %.2f fm^-3 along the QE peak nu = q^2/2M (non-relativistic)\n", rho0);
  std::printf("# q[GeV]  nu[GeV]  Q2[GeV2]   C_N    C_T    C_L\n");
  for (double q : {0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2}) {
    const double nu = q * q / (2 * M);
    const auto f = rpa.Evaluate(nu, q, rho0, 1.0);
    std::printf("%6.2f  %7.4f  %7.4f   %6.3f %6.3f %6.3f\n", q, nu, q * q - nu * nu, f.CN, f.CT, f.CL);
  }
}

void Momentum() {
  const Nucleus& C = Carbon();
  RelativisticFermiGas rfg;
  LocalFermiGas lfg;
  ShellModelSpectralFunction sf;
  Random rng(3);
  HitNucleonRequest req;
  req.species[0] = pdg::kNeutron;
  const int nb = 25;
  const double kmax = 0.5;
  const NuclearModel* models[3] = {&rfg, &lfg, &sf};
  const char* names[3] = {"RFG", "LFG", "SF"};
  std::vector<double> h[3];
  const long N = 2000000;
  std::printf("# 12C neutron momentum distribution k^2 n(k) (normalised to 1) [GeV^-1]\n");
  for (int m = 0; m < 3; ++m) {
    h[m].assign(nb + 1, 0.0);
    double sumER = 0, sumER2 = 0, sumK = 0, sumK2 = 0;
    for (long i = 0; i < N; ++i) {
      InitialState a;
      models[m]->Sample(C, req, a, rng);
      const double k = a.hit[0].p.Mag();
      h[m][std::min(static_cast<int>(k / kmax * nb), nb)] += 1;
      sumER += a.hit[0].removalEnergy;
      sumER2 += a.hit[0].removalEnergy * a.hit[0].removalEnergy;
      sumK += k;
      sumK2 += k * k;
    }
    const double mean = sumER / N;
    std::printf("# %s: <k> = %.4f GeV, rms k = %.4f GeV, <E_R> = %.2f MeV, rms(E_R) = %.2f MeV, k > %.2f: %.3f\n", names[m],
                sumK / N, std::sqrt(sumK2 / N), 1e3 * mean, 1e3 * std::sqrt(sumER2 / N - mean * mean), kmax, h[m][nb] / N);
  }
  std::printf("# k_lo  k_hi    RFG      LFG       SF\n");
  for (int b = 0; b < nb; ++b)
    std::printf("%5.3f %5.3f  %7.3f  %7.3f  %7.3f\n", b * kmax / nb, (b + 1) * kmax / nb, h[0][b] / N / (kmax / nb),
                h[1][b] / N / (kmax / nb), h[2][b] / N / (kmax / nb));
  // Removal-energy distribution of the SF (neutron), 2 MeV bins.
  std::vector<double> he(60, 0.0);
  for (long i = 0; i < N; ++i) {
    InitialState a;
    sf.Sample(C, req, a, rng);
    const int b = static_cast<int>(a.hit[0].removalEnergy / 0.002);
    if (b < 60) he[b] += 1;
  }
  std::printf("# SF neutron removal-energy distribution [per MeV]\n# E_lo[MeV] E_hi  P\n");
  for (int b = 0; b < 60; ++b) std::printf("#E %4d %4d  %.4f\n", 2 * b, 2 * b + 2, he[b] / N / 2);
  const NuclearDensity& d = NuclearDensity::For(C);
  std::printf("# density profile: r [fm], rho [fm^-3], k_F [GeV], V_n = T_F + S_n [MeV], gap Q - Q_LFG [MeV], V_out [MeV]\n");
  for (double r = 0; r <= 4.01; r += 0.5) {
    const double V = LocalFermiGas::Well(C, d, pdg::kNeutron, r);
    const double g = LocalFermiGas::EnergyGap(C, d, pdg::kNeutron, pdg::kProton, r);
    std::printf("#  %.1f  %.4f  %.4f  %6.2f  %6.2f  %6.2f\n", r, d.Rho(r), d.KFermi(pdg::kNeutron, r), 1e3 * V, 1e3 * g,
                1e3 * (V - g));
  }
}

struct Setup {
  const char* name;
  std::shared_ptr<const NuclearModel> model;
  std::shared_ptr<const QuasiElasticKinematics> kin;
};

std::vector<Setup> Setups() {
  QuasiElasticKinematics::Options o, orpa;
  orpa.rpa = true;
  return {{"RFG", std::make_shared<RelativisticFermiGas>(), std::make_shared<QuasiElasticKinematics>(FF(), o)},
          {"LFG", std::make_shared<LocalFermiGas>(), std::make_shared<QuasiElasticKinematics>(FF(), o)},
          {"LFG+RPA", std::make_shared<LocalFermiGas>(), std::make_shared<QuasiElasticKinematics>(FF(), orpa)},
          {"SF", std::make_shared<ShellModelSpectralFunction>(), std::make_shared<QuasiElasticKinematics>(FF(), o)},
          {"SF+RPA", std::make_shared<ShellModelSpectralFunction>(), std::make_shared<QuasiElasticKinematics>(FF(), orpa)}};
}

void Q2(double E, long N, int nu) {
  const Nucleus& C = Carbon();
  const auto setups = Setups();
  const std::vector<double> edges = {0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0};
  const int nb = static_cast<int>(edges.size()) - 1;
  std::vector<std::vector<double>> h(setups.size(), std::vector<double>(nb, 0));
  std::vector<double> tot(setups.size(), 0);
  HitNucleonRequest req;
  req.species[0] = nu > 0 ? pdg::kNeutron : pdg::kProton;
  for (size_t s = 0; s < setups.size(); ++s) {
    Random rng(7);
    Event ev;
    for (long i = 0; i < N; ++i) {
      InitialState init;
      setups[s].model->Sample(C, req, init, rng);
      const double w = setups[s].kin->Generate(nu, E, C, init, *setups[s].model, ev, rng);
      if (w <= 0) continue;
      tot[s] += w;
      for (int b = 0; b < nb; ++b)
        if (ev.Q2 >= edges[b] && ev.Q2 < edges[b + 1]) h[s][b] += w;
    }
  }
  std::printf("# CCQE %s on 12C at E = %.2f GeV: dsigma/dQ2 per nucleon [1e-38 cm^2/GeV^2]\n", nu > 0 ? "nu_mu" : "nubar_mu", E);
  std::printf("# sigma [1e-38 cm^2]:");
  for (size_t s = 0; s < setups.size(); ++s) std::printf("  %s %.4f", setups[s].name, tot[s] / N);
  std::printf("\n# Q2_lo Q2_hi  ");
  for (const auto& s : setups) std::printf(" %8s", s.name);
  std::printf("   LFG/RFG RPA/LFG  SF/RFG  RPA/SF\n");
  auto ratio = [](double a, double b) { return b > 0 ? a / b : 0.0; };
  for (int b = 0; b < nb; ++b) {
    const double w = (edges[b + 1] - edges[b]) * N;
    std::printf("%5.3f %5.3f  ", edges[b], edges[b + 1]);
    for (size_t s = 0; s < setups.size(); ++s) std::printf(" %8.4f", h[s][b] / w);
    std::printf("   %6.3f  %6.3f  %6.3f  %6.3f\n", ratio(h[1][b], h[0][b]), ratio(h[2][b], h[1][b]), ratio(h[3][b], h[0][b]),
                ratio(h[4][b], h[3][b]));
  }
}

void XSec(long N, int nu) {
  const Nucleus& C = Carbon();
  const auto setups = Setups();
  HitNucleonRequest req;
  req.species[0] = nu > 0 ? pdg::kNeutron : pdg::kProton;
  std::printf("# CCQE %s on 12C, sigma per %s [1e-38 cm^2]\n# E[GeV]", nu > 0 ? "nu_mu" : "nubar_mu", nu > 0 ? "neutron" : "proton");
  for (const auto& s : setups) std::printf(" %8s", s.name);
  std::printf("\n");
  for (double E : {0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0}) {
    std::printf("%6.2f", E);
    for (const auto& s : setups) {
      Random rng(11, static_cast<uint64_t>(E * 1e4));
      Event ev;
      double sum = 0;
      for (long i = 0; i < N; ++i) {
        InitialState init;
        s.model->Sample(C, req, init, rng);
        sum += s.kin->Generate(nu, E, C, init, *s.model, ev, rng);
      }
      std::printf(" %8.4f", sum / N);
    }
    std::printf("\n");
  }
}

// 12C(e,e') quasi-elastic in the impulse approximation, vector current only (F_A = F_P = 0):
//   d2sigma/dOmega/domega = sum_N < |M|^2 k' / (64 pi^2 E_e E_on E') * E'/(k q) * (1/2) >,
//   |M|^2 = (e^2/Q2)^2 * L.H / 2 with L.H = ContractQE(F_A = 0) / 4  (electron spin average included),
// averaged over the nuclear model's (r, k) and the azimuth of k around q, with the energy-conserving
// delta function fixing the polar angle of k: E'(k + q) = E_vtx(k) + omega.
void EQE(double E, double thetaDeg, long N) {
  const Nucleus& C = Carbon();
  RelativisticFermiGas rfg;
  LocalFermiGas lfg;
  ShellModelSpectralFunction sf;
  const NuclearModel* models[3] = {&rfg, &lfg, &sf};
  const NucleonFormFactors ff = FF();
  const double th = thetaDeg * kPi / 180;
  const double alpha = 1 / 137.035999;
  const double e2 = 4 * kPi * alpha;
  const double kGeV2toNb = 0.389379e6;  // GeV^-2 -> nb
  std::printf("# 12C(e,e') QE, E = %.3f GeV, theta = %.2f deg: d2sigma/dOmega/domega [nb/sr/GeV], vector current only\n",
              E, thetaDeg);
  std::printf("# omega[GeV]  q[GeV]      RFG        LFG         SF\n");
  for (double omega = 0.01; omega < std::min(E, 1.2); omega += 0.01) {
    const double Ep = E - omega;
    if (Ep <= 0) break;
    const FourVector k(E, 0, 0, E);
    const FourVector kp(Ep, Ep * std::sin(th), 0, Ep * std::cos(th));
    const FourVector qv = k - kp;
    const double q = qv.P(), Q2 = -qv.M2();
    const ThreeVector qhat = qv.Vect().Unit();
    ThreeVector e1 = qhat.Cross({0, 1, 0}).Unit(), e2v = qhat.Cross(e1);
    double res[3] = {0, 0, 0};
    for (int m = 0; m < 3; ++m) {
      const NuclearModel& model = *models[m];
      Random rng(13, static_cast<uint64_t>(omega * 1e5));
      double sum = 0;
      for (long i = 0; i < N; ++i) {
        for (int species : {pdg::kProton, pdg::kNeutron}) {
          HitNucleonRequest req;
          req.species[0] = species;
          InitialState init;
          model.Sample(C, req, init, rng);
          const double M = pdg::NucleonMass(species);
          const double kk = init.hit[0].p.Mag();
          const double eOn = init.hit[0].OnShellEnergy(M);
          const double eVtx = model.VertexEnergy(C, init, 0, species);
          const double Eout = eVtx + omega;
          // cos(alpha) between k and q from E'^2 = M^2 + k^2 + q^2 + 2 k q cos(alpha)
          if (kk <= 0) continue;
          const double ca = (Eout * Eout - M * M - kk * kk - q * q) / (2 * kk * q);
          if (std::abs(ca) > 1) continue;
          const double sa = std::sqrt(1 - ca * ca), phi = 2 * kPi * rng.Uniform();
          const ThreeVector kvec = qhat * (kk * ca) + e1 * (kk * sa * std::cos(phi)) + e2v * (kk * sa * std::sin(phi));
          InitialState moved = init;
          moved.hit[0].p = kvec;
          const ThreeVector pout = kvec + qv.Vect();
          if (model.PauliBlocked(C, species, pout, moved)) continue;
          const FourVector pOn(eOn, kvec), pOut(Eout, pout);
          const double Qt2 = -(pOut - pOn).M2();
          if (Qt2 <= 0) continue;
          // EM form factors of this nucleon from the Sachs form factors.
          const SachsFormFactors s = ff.Vector().Sachs(Qt2);
          const double tau = Qt2 / (4 * M * M);
          const double GE = species == pdg::kProton ? s.GEp : s.GEn, GM = species == pdg::kProton ? s.GMp : s.GMn;
          WeakFormFactors w;
          w.F1 = (GE + tau * GM) / (1 + tau);
          w.F2 = (GM - GE) / (1 + tau);
          const double LH = ContractQE(k, kp, pOn, pOut, w, M, +1) / 4;
          const double M2amp = e2 * e2 / (Q2 * Q2) * LH / 2;
          sum += M2amp * Ep / (64 * kPi * kPi * E * eOn * Eout) * Eout / (kk * q) * 0.5 *
                 (species == pdg::kProton ? C.Z : C.N());
        }
      }
      res[m] = sum / N * kGeV2toNb;
    }
    std::printf("%7.3f  %7.4f  %9.3f  %9.3f  %9.3f\n", omega, q, res[0], res[1], res[2]);
  }
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: nuflow-lfg-validate factors|momentum|q2 E [N] [nu]|xsec [N] [nu]|eqe E theta [N]\n");
    return 1;
  }
  const std::string mode = argv[1];
  if (mode == "factors") Factors();
  else if (mode == "momentum") Momentum();
  else if (mode == "q2") Q2(std::atof(argv[2]), argc > 3 ? std::atol(argv[3]) : 2000000, argc > 4 ? std::atoi(argv[4]) : 14);
  else if (mode == "xsec") XSec(argc > 2 ? std::atol(argv[2]) : 400000, argc > 3 ? std::atoi(argv[3]) : 14);
  else if (mode == "eqe") EQE(std::atof(argv[2]), std::atof(argv[3]), argc > 4 ? std::atol(argv[4]) : 200000);
  else {
    std::fprintf(stderr, "unknown mode %s\n", mode.c_str());
    return 1;
  }
  return 0;
}
