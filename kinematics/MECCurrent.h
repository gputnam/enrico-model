#pragma once
// Relativistic charged-current two-body meson-exchange current (MEC) for 2p2h excitations, from
//   I. Ruiz Simo, J.E. Amaro, M.B. Barbaro, A. De Pace, J.A. Caballero, T.W. Donnelly,
//   "Relativistic model of 2p-2h meson exchange currents in (anti)neutrino scattering",
//   J. Phys. G 44, 065105 (2017) [arXiv:1604.08423]  (referred to as RS17 below).
// The current is the weak pion-production amplitude of Hernandez, Nieves & Valverde (HNV, PRD 76,
// 033005 (2007)) without the nucleon-pole terms, with the pion absorbed on a second nucleon:
// seagull (contact), pion-in-flight, pion-pole and Delta(1232) forward/backward terms, vector and
// axial (RS17 Eqs. 13-48). Isospin is handled with the decomposition of RS17 Eq. (50),
//   j^mu = tau_+(1) J1^mu + tau_+(2) J2^mu + (I_V)_+ J3^mu,
// and the charge-channel effective currents of RS17 Table I. See theory/mec.md.

#include <complex>
#include <string>

#include "core/Dirac.h"
#include "core/FourVector.h"
#include "core/Random.h"
#include "kinematics/FormFactors.h"

namespace nuflow {

struct MECParameters {
  // Couplings (RS17 Sec. III): f_piNN/m_pi = g_A/(2 f_pi) (Goldberger-Treiman, RS17 Eq. 19), f* = 2.13.
  double gA = 1.26;
  double fPi = 0.093;     // [GeV]
  double fStar = 2.13;    // pi N Delta coupling
  double mRho = 0.7758;   // F_rho(k^2) = 1 / (1 - k^2/m_rho^2)  [GeV]
  double MDelta = 1.232;  // [GeV]
  // N-Delta form factors (HNV): C3V(Q2) = C3V0 / (1+Q2/MV^2)^2 / (1+Q2/(4 MV^2)),
  // C5A(Q2) = C5A0 / (1+Q2/MA^2)^2.  Only C3V and C5A are kept, as in RS17 Sec. V.
  double C3V0 = 2.13;
  double MV = 0.84;        // [GeV]
  double C5A0 = 1.00;      // fit IV of arXiv:1001.4416 (same N-Delta axial vertex as CCRES)
  double MADelta = 0.93;   // [GeV]
  // Delta propagator: RS17 keep only the real part of 1/(P^2 - M_D^2 + i M_D Gamma(P^2)) (Sec. V),
  // so that the Delta -> N pi decay (real pion emission, CCRES) is not double counted.
  bool realDeltaPropagator = true;
  // Monopole strong form factors (Lambda^2 - m_pi^2)/(Lambda^2 - k^2) at the pi NN and pi N Delta
  // vertices of the exchanged pion, placed as in De Pace et al., NPA 726 (2003) 303, Eqs. (16)-(19)
  // (seagull F_piNN^2, pion-in-flight F_piNN(k1) F_piNN(k2), Delta F_piNN F_piNDelta), with their
  // cut-offs 1.3 / 1.15 GeV. These are needed to reproduce the published RS17 responses (theory/mec.md).
  // 0 = off.
  double lambdaPiNN = 1.3;      // [GeV]
  double lambdaPiNDelta = 1.15; // [GeV]
  double axialScale = 1.0;   // multiplies every axial term (1 = physical; 0 for tests)
  double vectorScale = 1.0;  // multiplies every vector term (1 = physical; 0 for tests)
  double backgroundScale = 1.0;  // multiplies the seagull, pion-in-flight and pion-pole terms (tests)
};

// Charge channel of the 2p2h transition (initial pair -> final pair).
enum class MECChannel {
  kNuNPtoPP,     // nu:    n(1) p(2) -> p p           J_pp = J1 + J3, final-state exchange
  kNuNNtoNP,     // nu:    n n -> n(1') p(2')         J_np = J2 + J3, initial-state exchange
  kNubarPPtoNP,  // nubar: p p -> n(1') p(2')         J    = J1 - J3, initial-state exchange
  kNubarNPtoNN,  // nubar: n(1) p(2) -> n n           J    = J2 - J3, final-state exchange
};

class MECCurrent {
 public:
  using cplx = std::complex<double>;
  // Spin index s = s1' + 2 s2' + 4 s1 + 8 s2 (s = 0: spin up along z, 1: down).
  struct Amplitudes {
    cplx a[16][4];
  };
  // The three isospin-stripped currents of RS17 Eq. (50) for one labelling of the nucleon lines.
  struct Currents {
    cplx J1[16][4], J2[16][4], J3[16][4];
  };

  MECCurrent(MECParameters par, NucleonFormFactors nucleonFF) : par_(par), nucleonFF_(std::move(nucleonFF)) {}

  const MECParameters& Parameters() const { return par_; }
  std::string Name() const;

  // J1, J2, J3 (upper Lorentz index, covariant spinors ubar u = 2M) for line 1: i1 -> f1 and line 2:
  // i2 -> f2, absorbing the W four-momentum Q = f1 + f2 - i1 - i2. All nucleons on shell (mass M).
  void ComputeCurrents(const FourVector& Q, const FourVector& f1, const FourVector& i1, const FourVector& f2,
                       const FourVector& i2, Currents& out) const;

  // Antisymmetrised amplitude A^mu for the charge channel (direct minus exchange). The elementary
  // hadron tensor of RS17 Eqs. (57)-(63) is r^{mu nu} = 2 sum_spins A^mu* A^nu.
  void Antisymmetrized(MECChannel ch, const FourVector& Q, const FourVector& h1, const FourVector& h2,
                       const FourVector& p1, const FourVector& p2, Amplitudes& out) const;

  // Delta width Gamma(s) [GeV] in the Delta rest frame (HNV), zero below the N pi threshold.
  double DeltaWidth(double s) const;

 private:
  MECParameters par_;
  NucleonFormFactors nucleonFF_;
};

// Lepton tensor L^{mu nu} = sum over lepton spins of l^mu l^nu*, with l^mu = ubar(k') gamma^mu (1-gamma5) u(k)
// (neutrino, helicity +1) or vbar(k) gamma^mu (1-gamma5) v(k') (antineutrino, helicity -1).
void LeptonTensor(const FourVector& k, const FourVector& kp, double mLep, int helicity, std::complex<double> L[4][4]);

// sum_{mu nu} L^{mu nu} A_mu A_nu^*  for one hadronic amplitude A^mu (upper index).
double ContractLepton(const std::complex<double> L[4][4], const std::complex<double> A[4]);

// Monte Carlo estimate of the 2p2h hadron tensor of symmetric nuclear matter normalised to Z protons
// (RS17 Eq. 9, V = 3 pi^2 Z / kF^3), q along z, for one charge channel:
//   W^{mu nu} = V/(2pi)^9 int d3p1' d3h1 d3h2 M^4/(E1 E2 E1' E2') delta(E1'+E2'-E1-E2-omegaEff) Theta r^{mu nu},
// r^{mu nu} = 2 sum_spins A^mu* A^nu (direct - exchange), or 4 sum_spins J^mu* J^nu (direct only).
// Responses: R^CC = W^00, R^CL = -Re(W^03+W^30)/2, R^LL = W^33, R^T = W^11+W^22, R^T' = Im(W^12-W^21)/2.
void MECHadronTensor(const MECCurrent& cur, MECChannel ch, double q, double omegaEff, double kF, double Z, long n,
                     Random& rng, std::complex<double> W[4][4], bool directOnly = false);

}  // namespace nuflow
