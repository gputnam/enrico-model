#pragma once
// Long-range (RPA) nuclear correlations of the CCQE (1p1h) response in the local Fermi gas, after
// J. Nieves, J.E. Amaro and M. Valverde, Phys. Rev. C 70, 055503 (2004) [nucl-th/0408005], Sec. III A
// and Appendix A 2 (theory/local_fermi_gas.md, Sec. 4).
//
// The 1p1h hadron tensor is split by form-factor structure (F1^2, F1 F2, F2^2, F_A^2, F_A F_P, F_P^2,
// (F1+F2) F_A) and, in the frame with q along z, by component block (00, 0z, zz, transverse). Each
// (structure, block) is multiplied by one of the RPA factors
//     C_N = 1 / |1 - c0 f'(rho) U_N|^2      (scalar-isovector channel, ph only)
//     C_L = 1 / |1 - U V_l|^2               (spin-isospin longitudinal: pi exchange + g')
//     C_T = 1 / |1 - U V_t|^2               (spin-isospin transverse: rho exchange + g')
// with U = U_N + U_Delta the ph + Delta-h Lindhard function of symmetric nuclear matter at the local
// density, evaluated at the in-medium energy and momentum transfer (q~0, |q|).

#include <complex>
#include <string>

#include "core/FourVector.h"
#include "kinematics/FormFactors.h"

namespace nuflow {

class ValenciaRPA {
 public:
  struct Parameters {
    double gPrime = 0.63;        // Landau-Migdal g' of the S = T = 1 channel (Nieves et al. 2004)
    double f2Over4Pi = 0.08;     // pi N N coupling f^2 / 4 pi
    double lambdaPi = 1.2;       // pi N N form-factor cutoff [GeV]
    double cRho = 2.0;           // C_rho
    double lambdaRho = 2.5;      // rho N N form-factor cutoff [GeV]
    double mRho = 0.770;         // [GeV]
    double fStarOverF = 2.13;    // f*_{pi N Delta} / f_{pi N N}
    double mDelta = 1.232;       // [GeV]
    // Free Delta width [GeV] in the forward Delta-h propagator (PDG Breit-Wigner width 117 MeV): the
    // zero-width static Delta-h bubble has a pole at q0 = sqrt(M_D^2 + q^2) - M, which the 1p1h phase
    // space reaches at |q| > ~0.6 GeV (theory/spectral_function.md, Sec. 5).
    double gammaDelta = 0.117;
    double c0 = 0.380;           // [GeV fm^3]
    double f0pIn = 0.33;         // f'_0 at the centre (interior) ...
    double f0pEx = 0.45;         // ... and outside (exterior), Eq. (33)-(34)
    bool delta = true;           // include the Delta-h Lindhard function
    // Spreading width of the 1p1h states [GeV]: U_N is evaluated at the complex energy q0 + i Gamma/2.
    // It damps the zero-sound pole that the undamped nuclear-matter response develops just above the
    // particle-hole continuum at |q| < 0.1 GeV (theory/spectral_function.md, Sec. 5). 0 = undamped.
    double gammaPH = 0.005;
    // Range of validity of the non-relativistic effective interaction: the corrections are applied in
    // full for |q| <= qFull and faded linearly to zero at |q| = qZero (theory/local_fermi_gas.md, Sec. 4.4).
    double qFull = 0.6;          // [GeV] upper end of the data the model was validated on
    double qZero = 1.0;          // [GeV] the static Delta-h bubble approaches the QE region above this
  };

  struct Factors {
    double CN = 1, CT = 1, CL = 1;
  };

  ValenciaRPA() = default;
  explicit ValenciaRPA(Parameters p) : par_(p) {}
  const Parameters& Params() const { return par_; }
  std::string Describe() const;

  // Non-relativistic Lindhard function of symmetric nuclear matter (spin-isospin degeneracy 4, direct +
  // crossed ph bubbles) [GeV^2], for energy transfer q0 >= 0 and momentum |q| [GeV], Fermi momentum kF.
  static std::complex<double> LindhardN(double q0, double q, double kF, double M);
  // The same at complex energy q0 + i gamma/2 (analytic continuation, gamma > 0).
  static std::complex<double> LindhardN(double q0, double q, double kF, double M, double gamma);
  // Delta-h Lindhard function (static Delta, constant width gammaDelta in the forward term) with the
  // N-Delta coupling ratio built in [GeV^2]; rho in GeV^3.
  std::complex<double> LindhardDelta(double q0, double q, double rho) const;
  // Effective spin-isospin interactions V_l, V_t [GeV^-2] at (q0, |q|).
  double Vl(double q0, double q) const;
  double Vt(double q0, double q) const;

  // RPA factors at energy transfer q0 and momentum |q| of the particle-hole pair, local density rho
  // [fm^-3] and rhoRatio = rho / rho(centre), including the validity fade in |q|. q0 is clamped into the
  // 1p1h support of the local Fermi sea, [max(0, q^2 - 2 q k_F), q^2 + 2 q k_F] / 2M (a nucleon of a
  // spectral function may lie outside the local Fermi sphere; theory/spectral_function.md, Sec. 5).
  Factors Evaluate(double q0, double q, double rhoFm3, double rhoRatio) const;

  // L_{mu nu} H^{mu nu}_RPA for the CCQE vertex (same arguments and normalization as ContractQE). With
  // all factors equal to one it reproduces ContractQE.
  static double Contract(const FourVector& k, const FourVector& kp, const FourVector& p, const FourVector& pp,
                         const WeakFormFactors& ff, double M, int helicity, const Factors& f);

 private:
  Parameters par_;
};

}  // namespace nuflow
