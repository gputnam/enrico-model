// nuflow-mec-response: 2p2h MEC nuclear response functions R^CC, R^CL, R^LL, R^T, R^T' [GeV^-1] of
// symmetric nuclear matter (RFG) normalised to a nucleus with Z protons, in the conventions of
// Ruiz Simo et al., J. Phys. G 44, 065105 (2017) [arXiv:1604.08423], Eqs. (3)-(9), (57)-(63):
//   W^{mu nu} = V/(2pi)^9 int d3p1' d3h1 d3h2 M^4/(E1 E2 E1' E2') delta(E1'+E2'-E1-E2-omega') Theta r^{mu nu},
//   V = 3 pi^2 Z / kF^3, omega' = omega - eps (separation energy of the pair).
// It also prints the 1p1h RFG responses (one-body CC current, neutrons -> protons) in the same
// normalisation, for the 2p2h/1p1h comparison of RS17 Fig. 8. Used for theory/mec.md validation.
//
// usage: nuflow-mec-response [--q GeV] [--wmin GeV] [--wmax GeV] [--nw N] [--kf GeV] [--eps GeV]
//                            [--z Z] [--n samples] [--nubar] [--direct] [--full-delta] [--lambda GeV (piNN; 0 = no strong FF)]
//                            [--vector-only] [--axial-only] [--no-delta] [--delta-only]
//                            [--c5a C5A0] [--mad MA] [--channels pp|np|all]
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include <yaml-cpp/yaml.h>

#include "core/Constants.h"
#include "core/Dirac.h"
#include "core/Random.h"
#include "kinematics/FormFactors.h"
#include "kinematics/MECCurrent.h"

using namespace nuflow;
using namespace nuflow::constants;
using cplx = std::complex<double>;

namespace {

FourVector OnShell(double M, const ThreeVector& p) { return {std::sqrt(p.Mag2() + M * M), p}; }

ThreeVector InSphere(Random& r, double kF) {
  const double k = kF * std::cbrt(r.UniformOpen());
  const double c = 2 * r.Uniform() - 1, s = std::sqrt(std::max(0.0, 1 - c * c)), ph = 2 * kPi * r.Uniform();
  return {k * s * std::cos(ph), k * s * std::sin(ph), k * c};
}

struct Resp {
  double CC = 0, CL = 0, LL = 0, T = 0, Tp = 0;
};

Resp FromTensor(const cplx W[4][4]) {
  Resp r;
  r.CC = W[0][0].real();
  r.CL = -0.5 * (W[0][3] + W[3][0]).real();
  r.LL = W[3][3].real();
  r.T = (W[1][1] + W[2][2]).real();
  r.Tp = (cplx(0, -0.5) * (W[1][2] - W[2][1])).real();
  return r;
}

}  // namespace

int main(int argc, char** argv) {
  double q = 0.8, wmin = 0.02, wmax = -1, kF = 0.228, eps = 0.040, Z = 6, lambda = 1.3;
  int nw = 20;
  long n = 200000;
  bool nubar = false, direct = false, fullDelta = false;
  double c5a = 1.00, mad = 0.93;
  std::string chans = "all";
  double axScale = 1, vecScale = 1;
  bool noDelta = false, deltaOnly = false;
  for (int i = 1; i < argc; ++i) {
    auto next = [&]() { return std::atof(argv[++i]); };
    if (!std::strcmp(argv[i], "--q")) q = next();
    else if (!std::strcmp(argv[i], "--wmin")) wmin = next();
    else if (!std::strcmp(argv[i], "--wmax")) wmax = next();
    else if (!std::strcmp(argv[i], "--nw")) nw = static_cast<int>(next());
    else if (!std::strcmp(argv[i], "--kf")) kF = next();
    else if (!std::strcmp(argv[i], "--eps")) eps = next();
    else if (!std::strcmp(argv[i], "--z")) Z = next();
    else if (!std::strcmp(argv[i], "--n")) n = static_cast<long>(next());
    else if (!std::strcmp(argv[i], "--lambda")) lambda = next();
    else if (!std::strcmp(argv[i], "--c5a")) c5a = next();
    else if (!std::strcmp(argv[i], "--mad")) mad = next();
    else if (!std::strcmp(argv[i], "--nubar")) nubar = true;
    else if (!std::strcmp(argv[i], "--direct")) direct = true;
    else if (!std::strcmp(argv[i], "--full-delta")) fullDelta = true;
    else if (!std::strcmp(argv[i], "--vector-only")) axScale = 0;
    else if (!std::strcmp(argv[i], "--axial-only")) vecScale = 0;
    else if (!std::strcmp(argv[i], "--no-delta")) noDelta = true;
    else if (!std::strcmp(argv[i], "--delta-only")) deltaOnly = true;
    else if (!std::strcmp(argv[i], "--channels")) chans = argv[++i];
    else {
      std::fprintf(stderr, "unknown option %s\n", argv[i]);
      return 1;
    }
  }
  if (wmax < 0) wmax = q;

  MECParameters par;
  par.realDeltaPropagator = !fullDelta;
  par.lambdaPiNN = lambda;
  par.lambdaPiNDelta = lambda > 0 ? 1.15 : 0;  // --lambda 0 switches both off
  par.C5A0 = c5a;
  par.MADelta = mad;
  par.axialScale = axScale;
  par.vectorScale = vecScale;
  if (noDelta) par.fStar = 0;
  par.backgroundScale = deltaOnly ? 0 : 1;
  const NucleonFormFactors nff = NucleonFormFactors::FromConfig(YAML::Load("{vector: dipole, axial: dipole, gA: 1.26, MA: 1.05}"));
  const MECCurrent cur(par, nff);
  const double M = kMassNucleon;
  const double V = 3 * kPi * kPi * Z / (kF * kF * kF);
  const double sphere = 4 * kPi / 3 * kF * kF * kF;
  const double pref = V * sphere * sphere / std::pow(2 * kPi, 9);

  std::printf("# 2p2h MEC responses [GeV^-1], q = %.3f GeV, kF = %.3f GeV, eps = %.3f GeV, Z = %g, %s, %s Delta, Lambda = %g\n",
              q, kF, eps, Z, nubar ? "nubar" : "nu", fullDelta ? "full" : "real-part", lambda);
  std::printf("# omega      R_CC      R_CL      R_LL       R_T      R_T'   |  pp_T   np_T  | 1p1h: R_CC  R_CL  R_LL  R_T  R_T'\n");
  for (int iw = 0; iw < nw; ++iw) {
    const double w = nw == 1 ? wmin : wmin + (wmax - wmin) * iw / (nw - 1);
    const double wp = w - eps;
    const FourVector Q(wp, 0, 0, q);
    cplx W[2][4][4] = {};
    Random rng(1234, iw);
    for (long i = 0; i < n; ++i) {
      const FourVector h1 = OnShell(M, InSphere(rng, kF)), h2 = OnShell(M, InSphere(rng, kF));
      const FourVector P = Q + h1 + h2;
      const double s = P.M2();
      if (P.e <= 0 || s <= 4 * M * M) continue;
      const double W2 = std::sqrt(s);
      const double pstar = std::sqrt(s / 4 - M * M);
      const double c = 2 * rng.Uniform() - 1, sn = std::sqrt(std::max(0.0, 1 - c * c)), ph = 2 * kPi * rng.Uniform();
      const FourVector pcm(W2 / 2, pstar * sn * std::cos(ph), pstar * sn * std::sin(ph), pstar * c);
      const ThreeVector b = P.BoostVector();
      const FourVector p1 = pcm.Boosted(b), p2 = P - p1;
      if (p1.P2() < kF * kF || p2.P2() < kF * kF) continue;
      const double jac = 4 * kPi * pstar / W2 / (16 * h1.e * h2.e);
      for (int ic = 0; ic < 2; ++ic) {
        if (chans == "pp" && ic == 1) continue;
        if (chans == "np" && ic == 0) continue;
        const MECChannel ch = nubar ? (ic == 0 ? MECChannel::kNubarNPtoNN : MECChannel::kNubarPPtoNP)
                                    : (ic == 0 ? MECChannel::kNuNPtoPP : MECChannel::kNuNNtoNP);
        MECCurrent::Amplitudes A;
        if (direct) {
          MECCurrent::Currents d;
          cur.ComputeCurrents(Q, p1, h1, p2, h2, d);
          for (int sidx = 0; sidx < 16; ++sidx)
            for (int mu = 0; mu < 4; ++mu) {
              cplx v = 0;
              switch (ch) {
                case MECChannel::kNuNPtoPP: v = d.J1[sidx][mu] + d.J3[sidx][mu]; break;
                case MECChannel::kNuNNtoNP: v = d.J2[sidx][mu] + d.J3[sidx][mu]; break;
                case MECChannel::kNubarPPtoNP: v = d.J1[sidx][mu] - d.J3[sidx][mu]; break;
                case MECChannel::kNubarNPtoNN: v = d.J2[sidx][mu] - d.J3[sidx][mu]; break;
              }
              A.a[sidx][mu] = v;
            }
        } else {
          cur.Antisymmetrized(ch, Q, h1, h2, p1, p2, A);
        }
        // direct only: 4 sum |J|^2 (RS17 Eq. 59 without the exchange term) = 2 * 2 sum |J|^2
        const double fac = direct ? 4.0 : 2.0;
        for (int sidx = 0; sidx < 16; ++sidx)
          for (int mu = 0; mu < 4; ++mu)
            for (int nu = 0; nu < 4; ++nu) W[ic][mu][nu] += fac * jac * std::conj(A.a[sidx][mu]) * A.a[sidx][nu];
      }
    }
    cplx Wt[4][4];
    Resp rc[2];
    for (int ic = 0; ic < 2; ++ic) {
      cplx Wc[4][4];
      for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu) Wc[mu][nu] = W[ic][mu][nu] * (pref / n);
      rc[ic] = FromTensor(Wc);
    }
    for (int mu = 0; mu < 4; ++mu)
      for (int nu = 0; nu < 4; ++nu) Wt[mu][nu] = (W[0][mu][nu] + W[1][mu][nu]) * (pref / n);
    const Resp r = FromTensor(Wt);

    // 1p1h RFG (neutron -> proton for nu): W = V/(2pi)^2 int_0^kF h dh/(4 E q) sum |j_cov|^2 at cos fixed by energy.
    Resp r1;
    {
      const dirac::Gammas& G = dirac::Gammas::Get();
      const int nh = 400;
      cplx W1[4][4] = {};
      for (int ih = 0; ih < nh; ++ih) {
        const double h = kF * (ih + 0.5) / nh;
        const double E = std::sqrt(h * h + M * M), Ep = E + wp;
        if (Ep * Ep - M * M < kF * kF) continue;
        const double c = (Ep * Ep - M * M - h * h - q * q) / (2 * h * q);
        if (std::abs(c) > 1) continue;
        const FourVector hv(E, h * std::sqrt(1 - c * c), 0, h * c);
        const FourVector pv = hv + Q;
        const WeakFormFactors ff = nff.Evaluate(-Q.M2());
        const dirac::Mat Qs = G.Slash(Q);
        dirac::Mat S[4];
        for (int mu = 0; mu < 4; ++mu) {
          S[mu] = G.g[mu] * ff.F1;
          S[mu].AddScaled(G.g[mu] * Qs - Qs * G.g[mu], -ff.F2 / (4 * M));
          S[mu].AddScaled(G.g5g[mu], nubar ? -ff.FA : -ff.FA);
        }
        // spinors
        auto U = [&](const FourVector& p, int sp) {
          std::array<cplx, 4> u{};
          const double nn = std::sqrt(p.e + M), inv = 1 / (p.e + M);
          const cplx I(0, 1);
          if (sp == 0) u = {nn, 0, nn * p.pz * inv, nn * (p.px + I * p.py) * inv};
          else u = {0, nn, nn * (p.px - I * p.py) * inv, -nn * p.pz * inv};
          return u;
        };
        for (int sp = 0; sp < 2; ++sp)
          for (int si = 0; si < 2; ++si) {
            auto uo = U(pv, sp), ui = U(hv, si);
            std::array<cplx, 4> ub = {std::conj(uo[0]), std::conj(uo[1]), -std::conj(uo[2]), -std::conj(uo[3])};
            cplx j[4];
            for (int mu = 0; mu < 4; ++mu) {
              cplx t = 0;
              for (int a = 0; a < 4; ++a)
                for (int bb = 0; bb < 4; ++bb) t += ub[a] * S[mu](a, bb) * ui[bb];
              j[mu] = t;
            }
            for (int mu = 0; mu < 4; ++mu)
              for (int nu = 0; nu < 4; ++nu)
                W1[mu][nu] += std::conj(j[mu]) * j[nu] * (h / (4 * E * q)) * (kF / nh);
          }
      }
      for (int mu = 0; mu < 4; ++mu)
        for (int nu = 0; nu < 4; ++nu) W1[mu][nu] *= V / std::pow(2 * kPi, 2);
      r1 = FromTensor(W1);
    }
    std::printf("%7.3f %9.4f %9.4f %9.4f %9.4f %9.4f | %6.3f %6.3f | %7.3f %7.3f %7.3f %7.3f %7.3f\n", w, r.CC, r.CL,
                r.LL, r.T, r.Tp, rc[0].T, rc[1].T, r1.CC, r1.CL, r1.LL, r1.T, r1.Tp);
    std::fflush(stdout);
  }
  return 0;
}
