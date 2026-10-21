#pragma once
// Hadronic current for charged-current single-pion production  W+ N -> N' pi  in the model of
// E. Hernandez, J. Nieves, M. Valverde, Phys. Rev. D 76, 033005 (2007) [hep-ph/0701149] (HNV):
// direct and crossed Delta(1232) pole terms plus the leading-order chiral (non-linear sigma model)
// background: nucleon pole, crossed nucleon pole, contact, pion pole and pion-in-flight terms.
// N-Delta axial form factor C5A from fit IV of Hernandez et al., Phys. Rev. D 81, 085046 (2010)
// [arXiv:1001.4416] (ANL + BNL deuterium data, background and deuteron effects included).
// See theory/pion_production.md for the equations, conventions and validation.

#include <string>

#include "core/Dirac.h"
#include "core/FourVector.h"
#include "kinematics/FormFactors.h"

namespace nuflow {

struct HNVParameters {
  // N-Delta axial form factor (dipole), fit IV of arXiv:1001.4416.
  double C5A0 = 1.00;
  double MADelta = 0.93;  // [GeV]
  // N-Delta vector form factors, HNV Eq. (46) (Lalakulich-Paschos-Piranishvili).
  double C3V0 = 2.13, C4V0 = -1.51, C5V0 = 0.48;
  double MV = 0.84;  // [GeV]
  // Delta and couplings, HNV Sec. II C.
  double MDelta = 1.232;  // [GeV]
  double fStar = 2.14;    // pi N Delta coupling f*
  double fPi = 0.0924;    // pion decay constant [GeV]
  double gA = 1.26;       // nucleon axial coupling in the chiral background
  double mRho = 0.7758;   // rho mass in F_rho(t) [GeV]
  bool background = true; // include the non-resonant chiral terms (and the crossed Delta)
  double axialScale = 1.0; // multiplies every axial-current term (1 = physical; 0 is used to test CVC)
};

// N-Delta transition form factors C_i^{V,A}(Q2) of HNV Eq. (40).
struct NDeltaFormFactors {
  double C3V = 0, C4V = 0, C5V = 0, C3A = 0, C4A = 0, C5A = 0, C6A = 0;
};

// The seven HNV amplitudes (cos theta_C stripped), each a Dirac matrix per Lorentz index mu.
struct PionCurrentPieces {
  dirac::Mat dp[4], cdp[4], np[4], cnp[4], ct[4], pp[4], pf[4];
};

// Isospin combinations of the pieces (HNV Eq. 51 and Eqs. 10-13):
//   PPiPlus  = <p pi+ | j_cc+ | p>,   NPiPlus = <n pi+ | j_cc+ | n>,
//   PPiZero  = <p pi0 | j_cc+ | n> = -(PPiPlus - NPiPlus)/sqrt(2).
// Antineutrino amplitudes follow from Eqs. (11)-(13): <n pi- |j_cc-| n> = PPiPlus,
// <p pi- |j_cc-| p> = NPiPlus, <n pi0 |j_cc-| p> = -PPiZero.
enum class PionIsospin { PPiPlus, NPiPlus, PPiZero };

class PionProductionCurrent {
 public:
  PionProductionCurrent(HNVParameters par, NucleonFormFactors nucleonFF)
      : par_(par), nucleonFF_(std::move(nucleonFF)) {}

  const HNVParameters& Parameters() const { return par_; }
  std::string Name() const;

  // Delta width in its rest frame, HNV Eq. (45) [GeV]; zero below the N pi threshold.
  double DeltaWidth(double s) const;

  NDeltaFormFactors DeltaFormFactors(double Q2) const;
  WeakFormFactors NucleonFF(double Q2) const { return nucleonFF_.Evaluate(Q2); }

  // Pieces for on-shell initial nucleon p, final nucleon pp, pion kpi and q = pp + kpi - p.
  void Pieces(const FourVector& p, const FourVector& pp, const FourVector& kpi, PionCurrentPieces& out) const;
  // Current S^mu for an isospin channel built from the pieces.
  void Combine(const PionCurrentPieces& pieces, PionIsospin iso, dirac::Mat S[4]) const;

  // L_{mu nu} H^{mu nu} (full spin sums, no averaging) for the lepton momenta k (neutrino) and kp
  // (charged lepton, mass mLep), helicity +1 (nu) / -1 (nubar), and the hadronic current S with
  // initial/final nucleon momenta p (mass mIn) and pp (mass mOut).
  static double Contract(const FourVector& k, const FourVector& kp, double mLep, int helicity, const FourVector& p,
                         double mIn, const FourVector& pp, double mOut, const dirac::Mat S[4]);

 private:
  HNVParameters par_;
  NucleonFormFactors nucleonFF_;
};

}  // namespace nuflow
