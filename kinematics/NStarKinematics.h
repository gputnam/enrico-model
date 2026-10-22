#pragma once
// Charged-current excitation of the isospin-1/2 resonances of the second resonance region,
//   nu n -> l- R+,  nubar p -> l+ R0,   R = N(1440) P11, N(1520) D13, N(1535) S11,
// in the formalism and with the form factors of O. Lalakulich, E.A. Paschos, G. Piranishvili,
// Phys. Rev. D 74, 014009 (2006) [hep-ph/0602210] (LPP): vector form factors from the JLab/MAID
// electroproduction helicity amplitudes (CVC, isovector part C^V = C^n - C^p), axial couplings from
// PCAC and the R -> N pi widths. The resonances are added incoherently (no interference with each
// other or with the Delta + background of CCRES), with an energy-dependent Breit-Wigner, and decay
// to N pi, N eta (S11) or N pi pi with energy-dependent partial widths. Each resonance is excited
// only for W below w_max (= CCRES w_max = CCDIS w_min), where the inclusive Bodek-Yang continuum takes
// over. The hadron tensor is evaluated numerically from Dirac traces (spin-1/2 and Rarita-Schwinger
// spin-3/2); tests/test_nstar.cc checks it against the closed-form LPP structure functions.
// See theory/nstar.md.

#include <string>
#include <vector>

#include "core/FourVector.h"
#include "kinematics/KinematicsGenerator.h"

namespace nuflow {

struct NStarResonance {
  std::string name;  // "P11_1440", "D13_1520", "S11_1535"
  int twoJ;          // 1 or 3
  int parity;        // +1 or -1
  double M, Gamma0;  // Breit-Wigner mass and total width at W = M [GeV] (LPP Table I)
  int lPi;           // orbital angular momentum of the N pi decay
  double xPi, xEta, xPiPi;  // branching fractions at W = M (sum 1)
};

// Transition form factors at Q2 (weak, W+ n -> R+, LPP conventions).
struct NStarFormFactors {
  // spin 1/2: g1V, g2V, g1A, g3A;  spin 3/2: C3V, C4V, C5V, C5A, C6A (C3A = C4A = 0, LPP)
  double g1V = 0, g2V = 0, g1A = 0, g3A = 0;
  double C3V = 0, C4V = 0, C5V = 0, C5A = 0, C6A = 0;
};

class NStarKinematics : public KinematicsGenerator {
 public:
  struct Options {
    double wMax = 1.6;  // upper limit of W [GeV] (CCRES w_max = CCDIS w_min)
    double MV = 0.84;   // vector mass of the LPP form factors [GeV]
    double MA = 1.05;   // axial mass of the LPP form factors [GeV]
    double proposalQ2Scale = 1.0;
    double vectorScale = 1.0;  // multiplies all vector form factors (1 = physical; tests only)
    double axialScale = 1.0;   // multiplies all axial form factors (1 = physical; tests only)
    std::vector<std::string> resonances = {"P11_1440", "D13_1520", "S11_1535"};
  };

  explicit NStarKinematics(Options opt);

  std::string Name() const override;
  double Generate(int nuPdg, double enu, const Nucleus& nucleus, const InitialState& initial,
                  const NuclearModel& nuclearModel, Event& event, Random& rng) const override;

  // --- Physics building blocks (public for tests and validation tools).
  const std::vector<NStarResonance>& Resonances() const { return res_; }
  static NStarResonance Lookup(const std::string& name);
  NStarFormFactors FormFactors(const NStarResonance& r, double Q2) const;
  // Partial widths at mass W [GeV]: N pi, N eta, N pi pi; total = sum.
  void PartialWidths(const NStarResonance& r, double W, double& gPi, double& gEta, double& gPiPi) const;
  // Relativistic Breit-Wigner spectral function A(W^2) = (1/pi) M Gamma(W) / ((W^2-M^2)^2 + M^2 Gamma(W)^2)
  // [GeV^-2], normalised to ~1 over W^2.
  double Spectral(const NStarResonance& r, double W) const;
  // L_{mu nu} H^{mu nu} (full spin sums, no averaging) for W+ N -> R with nucleon momentum p (mass mIn),
  // resonance momentum pR (q = pR - p), resonance mass parameter M_R in the vertex / spin projector.
  double ContractLH(const NStarResonance& r, const FourVector& k, const FourVector& kp, double mLep, int helicity,
                    const FourVector& p, double mIn, const FourVector& pR) const;

 private:
  Options opt_;
  std::vector<NStarResonance> res_;
};

}  // namespace nuflow
