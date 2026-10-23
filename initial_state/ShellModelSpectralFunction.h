#pragma once
// Shell-model spectral function of 12C with short-range-correlated (SRC) pairs
// (theory/spectral_function.md).
//
//   P(k, E) = Z_MF * sum_{a = 1s, 1p} (n_a / N) |phi_a(k)|^2 L_a(E)  +  (1 - Z_MF) * P_SRC(k, E)
//
// * Mean field: harmonic-oscillator 1s and 1p orbitals with the oscillator length b of the HO point
//   density used by the cascade (fsi/NuclearDensity, fixed by the measured charge radius), IPSM
//   occupancies (2 and Z-2 per species), removal-energy distributions L_a(E) from 12C(e,e'p): Lorentzians
//   with the energy-dependent spreading width of Brown & Rho, Gamma(E) = a (E - E_F)^2 / (b + (E - E_F)^2),
//   1s centroid 39 MeV, a = 24 MeV, b = 590 MeV^2 (Lapikas et al., Phys. Rev. C 61, 064325 (2000)); the
//   1p strength is concentrated in the three low-lying 11B states (centroid 16.6 MeV, same width law).
//   Neutron removal energies are shifted by S_n - S_p (Coulomb).
// * SRC: a fraction 1 - Z_MF = 0.20 of the nucleons (Subedi et al., Science 320, 1476 (2008)) has
//   momentum k > k_SRC = 0.25 GeV with n(k) ~ 1/k^4 (contact scaling, Hen et al. PRC 92, 045205 (2015));
//   its partner carries -k + K_cm, K_cm Gaussian with sigma = 0.15 GeV per component (Cohen et al.,
//   PRL 121, 092501 (2018)), and is a nucleon of the other species with probability 0.95 (np dominance,
//   Hen et al., Science 346, 614 (2014)). The partner is emitted on shell and the A-2 residual is left
//   in its ground state, which fixes the (off-shell) energy of the struck nucleon.
// * The interaction point is sampled from the orbital's HO radial density (total density for SRC
//   nucleons); Pauli blocking uses the local k_F(r) of fsi/NuclearDensity at that point.
// Energies follow the asymptotic convention: the outgoing nucleon is produced on shell with p' = k + q
// and the struck nucleon has E = M_A - sqrt(M*_{A-1}^2 + k^2), M*_{A-1} = M_A - M_N + E_R.

#include <vector>

#include "initial_state/NuclearModel.h"

namespace nuflow {

class ShellModelSpectralFunction : public NuclearModel {
 public:
  struct Parameters {
    double mfFraction = 0.80;    // Z_MF; 1 - Z_MF nucleons in SRC pairs
    double e1s = 0.039;          // 1s centroid [GeV] (proton removal)
    double e1p = 0.0166;         // 1p centroid [GeV] (strength-weighted 11B 3/2-, 1/2-, 3/2- levels)
    double widthA = 0.024;       // Brown-Rho spreading width: a [GeV]
    double widthB = 590e-6;      // ... and b [GeV^2]
    double eMax = 0.100;         // upper end of the mean-field removal-energy distributions [GeV]
    double kSRC = 0.25;          // lower edge of the SRC momentum tail [GeV]
    double kMax = 1.0;           // upper edge [GeV]
    double sigmaCM = 0.15;       // pair centre-of-mass momentum width per component [GeV]
    double npFraction = 0.95;    // probability that the partner has the other isospin
    bool emitPartner = true;     // emit the SRC partner (false: it stays in the remnant)
  };

  ShellModelSpectralFunction();
  explicit ShellModelSpectralFunction(Parameters p);

  std::string Name() const override;
  void Sample(const Nucleus& nucleus, const HitNucleonRequest& request, InitialState& out,
              Random& rng) const override;
  bool PauliBlocked(const Nucleus& nucleus, int pdg, const ThreeVector& p) const override;
  bool PauliBlocked(const Nucleus& nucleus, int pdg, const ThreeVector& p, const InitialState& initial) const override;
  MediumPoint Medium(const Nucleus& nucleus, const InitialState& initial) const override;

  const Parameters& Params() const { return par_; }
  // Sample a proton removal energy [GeV] from shell s (0 = 1s, 1 = 1p), relative to nothing (absolute).
  double SampleRemovalEnergy(int shell, double u) const;

 private:
  void BuildTables();
  Parameters par_;
  // Inverse-CDF tables of the 1s and 1p removal-energy distributions (proton, 12C).
  std::vector<double> eGrid_[2];
};

}  // namespace nuflow
