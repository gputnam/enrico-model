#include "initial_state/NuclearModel.h"

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
  out.remnant = rem;
  out.remnantPdg = A > 0 ? pdg::IonCode(Z, A) : 0;
}

}  // namespace nuflow
