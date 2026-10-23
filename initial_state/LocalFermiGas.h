#pragma once
// Local Fermi gas (LFG) for the primary vertex (theory/local_fermi_gas.md).
//
// The interaction point r is sampled from the nucleon density rho(r) of fsi/NuclearDensity (the same
// harmonic-oscillator / 2pF profile the intranuclear cascade uses), and the struck nucleon's momentum
// uniformly in the local Fermi sphere k < k_F^N(r) = [3 pi^2 rho_N(r)]^(1/3). Pauli blocking of the
// outgoing nucleon uses the local k_F of its species at the same r.
//
// Energy balance (Nieves, Amaro & Valverde, Phys. Rev. C 70, 055503 (2004), Eq. 43): in the medium both
// nucleons are on shell in a common local potential, and the 1p1h transition costs the energy gap
//     Delta(r) = Q - Q_LFG(r),   Q_LFG(r) = E_F^final(r) - E_F^hit(r),   E_F = sqrt(M^2 + k_F^2),
// with Q the ground-state Q value of the CC transition (12C -> 12N: 16.827 MeV), i.e.
//     E_final(k + q) = E_hit(k) + q0 - Delta(r).
// The hit nucleon is bound in a well of depth V_hit(r) = T_F^hit(r) + S_hit (S = separation energy), so
// its removal energy is E_R = V_hit(r) - T(k) in [S_hit, S_hit + T_F(r)]; the outgoing nucleon then
// leaves a well of depth V_out(r) = V_hit(r) - Delta(r) (12C, nu: T_F(r) + 0.60 MeV = T_F + S_p(12N)),
// decelerating along its direction, so that asymptotic energies satisfy q0 = T'_asym + E_R.

#include "initial_state/NuclearModel.h"

namespace nuflow {

class NuclearDensity;

class LocalFermiGas : public NuclearModel {
 public:
  enum class EnergyBalance {
    kInMedium,    // Nieves et al.: on-shell in-medium nucleons, Q-value gap, outgoing nucleon leaves the well
    kAsymptotic,  // removal energy E_R at the vertex, outgoing nucleon asymptotic (no outgoing potential)
  };
  LocalFermiGas() = default;
  explicit LocalFermiGas(EnergyBalance e) : energy_(e) {}
  std::string Name() const override {
    return energy_ == EnergyBalance::kInMedium ? "LFG(HO/2pF density, Q-value gap)" : "LFG(HO/2pF density, asymptotic E_R)";
  }
  EnergyBalance Energy() const { return energy_; }
  void Sample(const Nucleus& nucleus, const HitNucleonRequest& request, InitialState& out,
              Random& rng) const override;
  // Without an interaction point the LFG cannot decide; this overload throws std::logic_error.
  bool PauliBlocked(const Nucleus& nucleus, int pdg, const ThreeVector& p) const override;
  bool PauliBlocked(const Nucleus& nucleus, int pdg, const ThreeVector& p, const InitialState& initial) const override;
  double VertexEnergy(const Nucleus& nucleus, const InitialState& initial, int i, int finalPdg) const override;
  double OutgoingPotential(const Nucleus& nucleus, const InitialState& initial, int finalPdg) const override;
  MediumPoint Medium(const Nucleus& nucleus, const InitialState& initial) const override;

  // Depth of the well binding a nucleon of species pdg at radius r: T_F(r) + S [GeV].
  static double Well(const Nucleus& nucleus, const NuclearDensity& density, int pdg, double r);
  // Energy gap Delta(r) = Q - Q_LFG(r) of the transition hitPdg -> finalPdg [GeV] (0 if no charge change).
  static double EnergyGap(const Nucleus& nucleus, const NuclearDensity& density, int hitPdg, int finalPdg, double r);

 private:
  EnergyBalance energy_ = EnergyBalance::kInMedium;
};

}  // namespace nuflow
