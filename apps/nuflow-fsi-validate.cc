// nuflow-fsi-validate: hadron-nucleus validation of the intranuclear cascade (theory/fsi_cascade.md).
//
//   nuflow-fsi-validate xsec                       free NN and pi N cross sections used by the cascade
//   nuflow-fsi-validate pion [pdg] [N]             pi-12C reaction / absorption / charge exchange / inelastic
//   nuflow-fsi-validate transparency [N]           proton transparency of 12C vs proton kinetic energy
//   nuflow-fsi-validate events CONFIG [N]          fates of primary hadrons in neutrino events (CC1pi+ etc.)

#include <cmath>
#include <cstdio>
#include <map>
#include <string>

#include "core/Config.h"
#include "core/Constants.h"
#include "core/PDG.h"
#include "fsi/HadronNucleonXSec.h"
#include "fsi/IntranuclearCascade.h"
#include "orchestration/Generator.h"
#include "orchestration/OutputWriter.h"

using namespace nuflow;

namespace {

const Nucleus& Carbon() { return NuclearDatabase::Instance().Get(1000060120); }

int RunXSec() {
  std::printf("# NN (Cugnon et al. 1996): plab[GeV]  pp_tot pp_el  np_tot np_el [mb]  B_pp B_np [GeV^-2]\n");
  for (double p : {0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0})
    std::printf("%5.2f  %7.1f %7.1f  %7.1f %7.1f  %5.2f %5.2f\n", p, hadxs::NNTotal(true, p), hadxs::NNElastic(true, p),
                hadxs::NNTotal(false, p), hadxs::NNElastic(false, p), hadxs::NNSlope(true, p), hadxs::NNSlope(false, p));
  std::printf("# free pi N (nucleon at rest): Tpi[MeV]  W  Gamma_Delta  sigma(pi+p)  sigma(pi-p el)  sigma(pi-p cx) [mb]\n");
  const double M = constants::kMassProton, m = constants::kMassPionCharged;
  for (double T : {0.03, 0.05, 0.085, 0.125, 0.165, 0.19, 0.205, 0.245, 0.315, 0.4, 0.5, 0.7, 1.0}) {
    const double E = T + m, W = std::sqrt(M * M + m * m + 2 * M * E);
    const double q = hadxs::CMMomentum(W, M, m), G = hadxs::DeltaWidth(W), wD = hadxs::DeltaRegionWeight(W);
    const double bw = (0.25 * G * G) / ((W - hadxs::kMassDelta) * (W - hadxs::kMassDelta) + 0.25 * G * G);
    const double s33 = wD * 8 * constants::kPi / (q * q) * 0.389379 * bw;
    const double bg = (1 - wD) * hadxs::kPiNHighEnergy;
    std::printf("%6.0f  %6.3f  %6.3f  %7.1f  %7.1f  %7.1f\n", T * 1e3, W, G, s33 + bg, s33 / 9 + bg, 2 * s33 / 9);
  }
  std::printf("# Oset-Salcedo -Im Sigma_Delta at rho0 [MeV]: Tpi  C_Q C_A2 C_A3\n");
  for (double T : {0.05, 0.085, 0.125, 0.165, 0.205, 0.245, 0.315}) {
    const auto s = hadxs::DeltaSelfEnergy(T, 1.0);
    std::printf("%6.0f  %6.1f %6.1f %6.1f\n", T * 1e3, s.cQ * 1e3, s.cA2 * 1e3, s.cA3 * 1e3);
  }
  const NuclearDensity d(Carbon());
  std::printf("# 12C density: HO b = %.3f fm, point rms = %.3f fm, rho(0) = %.4f fm^-3, kF(0) = %.3f GeV, Rmax = %.2f fm\n",
              d.B(), d.PointRms(), d.Rho(0), d.KFermi(pdg::kProton, 0), d.RMax());
  for (double r = 0; r <= 5.01; r += 0.5)
    std::printf("  r=%.1f rho=%.4f kF=%.3f\n", r, d.Rho(r), d.KFermi(pdg::kProton, r));
  return 0;
}

int RunPion(int pdgPi, long N) {
  const Nucleus& C = Carbon();
  IntranuclearCascade inc;
  const auto& d = IntranuclearCascade::DensityFor(C);
  const double R = d.RMax();
  const double area = constants::kPi * R * R * 10.0;  // fm^2 -> mb
  std::printf("# pi(%d) + 12C, %ld pions per energy, impact parameters uniform on a disk of R = %.2f fm\n", pdgPi, N, R);
  std::printf("# Tpi[MeV] sigma_reac  sigma_abs  sigma_cx  sigma_inel  sigma_abs+cx [mb]   <n_coll|reac>\n");
  Random rng(77);
  for (double T : {0.05, 0.085, 0.125, 0.165, 0.205, 0.245, 0.315, 0.5, 1.0}) {
    long nReac = 0, nAbs = 0, nCX = 0, nInel = 0;
    double nColl = 0;
    for (long i = 0; i < N; ++i) {
      double bx, by;
      do {
        bx = R * (2 * rng.Uniform() - 1);
        by = R * (2 * rng.Uniform() - 1);
      } while (bx * bx + by * by > R * R);
      IntranuclearCascade::Hadron h;
      h.pdg = pdgPi;
      h.x = {bx, by, -R - 0.01};
      h.dir = {0, 0, 1};
      h.estar = T;
      std::vector<IntranuclearCascade::Hadron> esc;
      int z = C.Z, n = C.N();
      IntranuclearCascade::Stats st;
      inc.Cascade(C, d, {h}, esc, z, n, rng, &st);
      const int nInt = st.pionQE + st.pionCX + st.pionAbs2 + st.pionAbs3;
      if (nInt == 0) continue;
      ++nReac;
      nColl += nInt;
      int nSame = 0, nPi0 = 0, nPiCh = 0;
      for (const auto& e : esc) {
        if (e.pdg == pdgPi) ++nSame;
        if (e.pdg == pdg::kPiZero) ++nPi0;
        if (std::abs(e.pdg) == 211) ++nPiCh;
      }
      if (nPiCh == 0 && nPi0 == 0) ++nAbs;
      else if (nPiCh == 0 && nPi0 == 1) ++nCX;
      else if (nSame >= 1) ++nInel;
    }
    std::printf("%6.0f  %8.1f  %8.1f  %8.1f  %8.1f  %8.1f   %.2f\n", T * 1e3, area * nReac / N, area * nAbs / N,
                area * nCX / N, area * nInel / N, area * (nAbs + nCX) / N, nReac ? nColl / nReac : 0.0);
  }
  return 0;
}

int RunTransparency(long N) {
  const Nucleus& C = Carbon();
  IntranuclearCascade inc;
  const auto& d = IntranuclearCascade::DensityFor(C);
  Random rng(78);
  std::printf("# proton knocked out at a density-weighted point, isotropic direction; transparency = P(no unblocked NN collision)\n");
  std::printf("# Tp[GeV]  Q2=2M*Tp[GeV2]  T  P(escape as proton)  <n_out nucleons>\n");
  for (double T : {0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0}) {
    long nClean = 0, nProton = 0;
    double nOut = 0;
    for (long i = 0; i < N; ++i) {
      IntranuclearCascade::Hadron h;
      h.pdg = pdg::kProton;
      h.x = d.SamplePosition(rng);
      const double c = 2 * rng.Uniform() - 1, s = std::sqrt(1 - c * c), phi = 2 * constants::kPi * rng.Uniform();
      h.dir = {s * std::cos(phi), s * std::sin(phi), c};
      h.estar = T;
      std::vector<IntranuclearCascade::Hadron> esc;
      int z = C.Z - 1, n = C.N();
      IntranuclearCascade::Stats st;
      inc.Cascade(C, d, {h}, esc, z, n, rng, &st);
      if (st.nucleonCollisions - st.nucleonPauliBlocked == 0) ++nClean;
      for (const auto& e : esc)
        if (e.pdg == pdg::kProton || e.pdg == pdg::kNeutron) nOut += 1;
      bool hasP = false;
      for (const auto& e : esc) hasP |= e.pdg == pdg::kProton;
      nProton += hasP;
    }
    std::printf("%5.2f  %5.2f  %.3f  %.3f  %.2f\n", T, 2 * constants::kMassProton * T, double(nClean) / N,
                double(nProton) / N, nOut / N);
  }
  return 0;
}

struct Collect : OutputWriter {
  std::map<std::string, double> c;
  long n = 0;
  double maxDE = 0;
  int maxQ = 0, maxB = 0;
  void Begin(const RunInfo&) override {}
  void End(const RunSummary&) override {}
  void Write(const Event& ev) override {
    ++n;
    FourVector in, out;
    int qIn = 0, qOut = 0, bIn = 0, bOut = 0;
    int nPipPre = 0, nPiPre = 0, nPipFin = 0, nPi0Fin = 0, nPimFin = 0, nPFin = 0, nNFin = 0, nNucPre = 0;
    for (const auto& p : ev.particles) {
      if (p.status == kStatusBeam || p.status == kStatusTarget) {
        in = in + p.p;
        if (p.status == kStatusTarget) qIn += pdg::IonZ(p.pdg), bIn += pdg::IonA(p.pdg);
      }
      if (p.status == kStatusFinal || p.status == kStatusRemnantPostFSI) {
        out = out + p.p;
        if (p.status == kStatusRemnantPostFSI) {
          qOut += pdg::IonZ(p.pdg), bOut += pdg::IonA(p.pdg);
          continue;
        }
        if (p.pdg == pdg::kProton) ++qOut, ++bOut, ++nPFin;
        if (p.pdg == pdg::kNeutron) ++bOut, ++nNFin;
        if (p.pdg == 211) ++qOut, ++nPipFin;
        if (p.pdg == -211 || p.pdg == 13) --qOut;
        if (p.pdg == -211) ++nPimFin;
        if (p.pdg == -13) ++qOut;
        if (p.pdg == 111) ++nPi0Fin;
      }
      if (p.status == kStatusPreFSI) {
        nPipPre += p.pdg == 211;
        nPiPre += std::abs(p.pdg) == 211 || p.pdg == 111;
        nNucPre += p.pdg == pdg::kProton || p.pdg == pdg::kNeutron;
      }
    }
    maxDE = std::max({maxDE, std::abs(in.e - out.e), (in.Vect() - out.Vect()).Mag()});
    maxQ = std::max(maxQ, std::abs(qIn - qOut));
    maxB = std::max(maxB, std::abs(bIn - bOut));
    const std::string ch = std::to_string(ev.channel);
    c["n_" + ch] += 1;
    c["nP_" + ch] += nPFin;
    c["nN_" + ch] += nNFin;
    if (ev.channel == 400 && nPipPre == 1 && ev.targetPdg != pdg::kProton) {
      c["res_pip"] += 1;
      if (nPipFin + nPimFin + nPi0Fin == 0) c["res_pip_abs"] += 1;
      else if (nPipFin == 0 && nPi0Fin >= 1 && nPimFin == 0) c["res_pip_cx"] += 1;
      else if (nPipFin == 1 && nPi0Fin == 0 && nPimFin == 0) c["res_pip_1pip"] += 1;
      else c["res_pip_other"] += 1;
    }
    if (ev.targetPdg != pdg::kProton && nPiPre == 0 && nPipFin + nPimFin + nPi0Fin > 0) c["pion_from_nn_" + ch] += 1;
  }
};

int RunEvents(const std::string& cfg, long N) {
  YAML::Node config = LoadConfig(cfg);
  Generator gen(config);
  Collect col;
  gen.Run(N, {&col});
  std::printf("# %ld events of %s\n", col.n, cfg.c_str());
  std::printf("max |Delta E|, |Delta p| = %.3g GeV, max |Delta Q| = %d, max |Delta B| = %d\n", col.maxDE, col.maxQ,
              col.maxB);
  for (const char* ch : {"200", "300", "400", "600"}) {
    const double n = col.c["n_" + std::string(ch)];
    if (n <= 0) continue;
    std::printf("channel %s: %.0f events, <n_p> = %.3f, <n_n> = %.3f, pion made in FSI (no primary pion): %.4f\n", ch, n,
                col.c["nP_" + std::string(ch)] / n, col.c["nN_" + std::string(ch)] / n,
                col.c["pion_from_nn_" + std::string(ch)] / n);
  }
  const double r = col.c["res_pip"];
  if (r > 0)
    std::printf("CCRES on C with a primary pi+ (%.0f): absorbed %.3f, charge exchanged (pi0, no pi+-) %.3f, "
                "single pi+ %.3f, other %.3f\n", r, col.c["res_pip_abs"] / r, col.c["res_pip_cx"] / r,
                col.c["res_pip_1pip"] / r, col.c["res_pip_other"] / r);
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  const std::string mode = argc > 1 ? argv[1] : "";
  try {
    if (mode == "xsec") return RunXSec();
    if (mode == "pion") return RunPion(argc > 2 ? std::stoi(argv[2]) : 211, argc > 3 ? std::stol(argv[3]) : 200000);
    if (mode == "transparency") return RunTransparency(argc > 2 ? std::stol(argv[2]) : 200000);
    if (mode == "events" && argc > 2) return RunEvents(argv[2], argc > 3 ? std::stol(argv[3]) : 200000);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  std::fprintf(stderr, "usage: nuflow-fsi-validate xsec | pion [pdg] [N] | transparency [N] | events CONFIG [N]\n");
  return 2;
}
