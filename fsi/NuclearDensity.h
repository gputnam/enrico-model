#pragma once
// Nucleon number density of a nucleus, used by the intranuclear cascade (theory/fsi_cascade.md).
//
//  * A <= 16 (s and p shell, N = Z): harmonic-oscillator shell-model point density
//        rho(r) = rho0 (1 + c x^2) exp(-x^2),  x = r/b,  c = (Z-2)/3,
//    with b fixed by the measured rms charge radius after removing the nucleon finite size
//    (<r^2>_pt = r_ch^2 - r_p^2 - (N/Z) r_n^2 - 3/(4M^2)).
//  * heavier nuclei: two-parameter Fermi distribution with diffuseness a = 0.52 fm and the
//    half-density radius fixed by the same rms constraint.
// Protons and neutrons follow the same profile, scaled by Z/A and N/A.

#include "core/FourVector.h"
#include "core/Nucleus.h"
#include "core/Random.h"

namespace nuflow {

class NuclearDensity {
 public:
  explicit NuclearDensity(const Nucleus& nucleus);

  double Rho(double r) const;                       // total nucleon density [fm^-3]
  double RhoSpecies(int nucleonPdg, double r) const;  // proton or neutron density [fm^-3]
  // Local Fermi momentum of a species, k_F = (3 pi^2 rho_species)^(1/3) [GeV].
  double KFermi(int nucleonPdg, double r) const;
  double RMax() const { return rMax_; }             // radius where rho = 1e-4 rho(0) [fm]
  double RhoCentral() const { return Rho(0); }
  ThreeVector SamplePosition(Random& rng) const;    // r^2 rho(r) dr dOmega
  bool HarmonicOscillator() const { return ho_; }
  double B() const { return b_; }                   // HO length or 2pF radius [fm]
  double PointRms() const { return rmsPoint_; }     // rms radius of the point-nucleon density [fm]

 private:
  double Shape(double r) const;  // un-normalised profile
  int Z_ = 0, A_ = 0;
  bool ho_ = true;
  double b_ = 1, c_ = 0, a_ = 0.52, norm_ = 1, rMax_ = 5, rmsPoint_ = 0, maxR2Rho_ = 0;
};

}  // namespace nuflow
