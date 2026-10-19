#pragma once
// Free-nucleon charged-current quasi-elastic matrix element.
//
// ContractQE returns the spin-summed contraction  L_{mu nu} H^{mu nu}  of
//   L^{mu nu} = Tr[(k'slash + m) gamma^mu (1 - gamma5) kslash gamma^nu (1 - gamma5)]   (neutrino)
//   H^{mu nu} = Tr[(p'slash + M) Gamma^mu (pslash + M) Gammabar^nu]
// for on-shell nucleon momenta p, p' of common mass M and form factors evaluated at
// Q2 = -(p' - p)^2. The closed form is derived in theory/quasielastic.md (Sec. 3).
// The spin-averaged squared amplitude is  <|M|^2> = G_F^2 |V_ud|^2 / 4 * ContractQE(...).

#include "core/FourVector.h"
#include "kinematics/FormFactors.h"

namespace nuflow {

// helicity = +1 for neutrinos, -1 for antineutrinos.
double ContractQE(const FourVector& k, const FourVector& kp, const FourVector& p, const FourVector& pp,
                  const WeakFormFactors& ff, double M, int helicity);

// Llewellyn Smith (Phys. Rep. 3 (1972) 261) closed form dsigma/dQ2 [GeV^-4] for a free nucleon at rest
// (mass M), as restated e.g. in Formaggio & Zeller, Rev. Mod. Phys. 84, 1307 (2012), Eqs. (57)-(60).
double LlewellynSmithDSigmaDQ2(double enu, double Q2, double mLepton, double M, const WeakFormFactors& ff,
                               int helicity);

}  // namespace nuflow
