#include "initial_state/NuclearModel.h"

#include <cmath>

#include "core/Constants.h"
#include "core/PDG.h"

namespace nuflow {

void FillRemnant(const Nucleus& nucleus, InitialState& out) {
  // Remnant = target at rest minus the hit nucleons. It carries the removal energy as excitation,
  // so four-momentum is conserved exactly at the primary vertex.
  FourVector rem(nucleus.mass, 0, 0, 0);
  int Z = nucleus.Z, A = nucleus.A;
  for (int i = 0; i < out.nHit; ++i) {
    rem = rem - out.hit[i].Momentum();
    A -= 1;
    if (out.hit[i].pdg == pdg::kProton) Z -= 1;
  }
  if (out.hasSpectator) {
    rem = rem - out.spectator.Momentum();
    A -= 1;
    if (out.spectator.pdg == pdg::kProton) Z -= 1;
  }
  out.remnant = rem;
  out.remnantPdg = A > 0 ? pdg::IonCode(Z, A) : 0;
}

MediumPoint NuclearModel::Medium(const Nucleus& nucleus, const InitialState&) const {
  // rho = (kF_p^3 + kF_n^3) / (3 pi^2), converted to fm^-3.
  constexpr double kHbarC = 0.1973269804;
  const double kp = nucleus.kFermiProton / kHbarC, kn = nucleus.kFermiNeutron / kHbarC;
  MediumPoint m;
  m.rho = (kp * kp * kp + kn * kn * kn) / (3 * constants::kPi * constants::kPi);
  m.rhoCentral = m.rho;
  return m;
}

}  // namespace nuflow
