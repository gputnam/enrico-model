#pragma once
// Free hadron-nucleon cross sections and in-medium Delta parameters used by the intranuclear
// cascade (theory/fsi_cascade.md). Units: mb for cross sections, GeV otherwise.
//
//  * NN: parametrisations of the measured pp and np total and elastic cross sections and of the
//    elastic diffraction slope, J. Cugnon, D. L'Hote, J. Vandermeulen, Nucl. Instrum. Meth. B 111,
//    215 (1996) (functions of the beam momentum p_lab in the rest frame of the target nucleon).
//  * piN: Delta(1232) Breit-Wigner with the energy-dependent P-wave width of the CCRES channel
//    (f* = 2.14, M_Delta = 1.232 GeV), isospin Clebsch-Gordan factors, and above W ~ 1.4 GeV an
//    isospin-averaged non-resonant total cross section of 30 mb (PDG piN compilation).
//  * In-medium Delta: imaginary part of the Delta self-energy of E. Oset and L.L. Salcedo,
//    Nucl. Phys. A 468, 631 (1987): -Im Sigma = C_Q (rho/rho0)^alpha + C_A2 (rho/rho0)^beta
//    + C_A3 (rho/rho0)^(2 beta), rho0 = 0.17 fm^-3, coefficients quadratic in x = T_pi/m_pi
//    (valid for 85 < T_pi < 315 MeV; frozen at the end points outside).

namespace nuflow::hadxs {

// --- nucleon-nucleon (like = pp or nn, otherwise np) ---
double NNTotal(bool like, double plab);
double NNElastic(bool like, double plab);
double NNSlope(bool like, double plab);  // B in dsigma/dt ~ exp(B t) [GeV^-2]

// --- pion-nucleon ---
inline constexpr double kMassDelta = 1.232;
inline constexpr double kFStar = 2.14;
double DeltaWidth(double W);              // free Pi N width Gamma(W) [GeV]
double CMMomentum(double W, double m1, double m2);
// Fraction of the Delta-dominated description in use at W (1 below 1.35 GeV, 0 above 1.45 GeV).
double DeltaRegionWeight(double W);
inline constexpr double kPiNHighEnergy = 30.0;  // mb, isospin-averaged sigma_tot(pi N) above the Delta
inline constexpr double kPiNHighEnergySlope = 7.0;  // GeV^-2

// Squared isospin Clebsch-Gordan coefficient for pi(charge qpi) N(charge qN) -> Delta.
double DeltaIsospinFactor(int qpi, int qN);

struct OsetSelfEnergy {
  double cQ = 0, cA2 = 0, cA3 = 0;  // -Im Sigma_Delta pieces [GeV]
};
OsetSelfEnergy DeltaSelfEnergy(double tpi, double rhoOverRho0);
inline constexpr double kRho0 = 0.17;  // fm^-3

}  // namespace nuflow::hadxs
