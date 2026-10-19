#pragma once
// Global relativistic Fermi gas (Smith & Moniz, Nucl. Phys. B43 (1972) 605), with an optional
// high-momentum tail n(k) ~ 1/k^4 above k_F. The Bodek-Ritchie tail (Phys. Rev. D 23, 1070 (1981),
// Eq. 6) puts a fraction 6 (k_F a / pi)^2, a = 2 GeV^-1, of the nucleons on [k_F, 4 GeV].

#include "initial_state/NuclearModel.h"

namespace nuflow {

class RelativisticFermiGas : public NuclearModel {
 public:
  struct Options {
    bool bodekRitchie = false;  // Bodek-Ritchie tail (overrides tailFraction / tailCutoff)
    double tailFraction = 0.0;  // fraction of nucleons in the k > k_F tail (0 = pure Smith-Moniz)
    double tailCutoff = 0.5;    // maximum tail momentum [GeV]
  };
  static constexpr double kBodekRitchieA = 2.0;       // correlation range parameter a [GeV^-1]
  static constexpr double kBodekRitchieCutoff = 4.0;  // [GeV]

  RelativisticFermiGas() = default;
  explicit RelativisticFermiGas(Options opt) : opt_(opt) {}

  std::string Name() const override;
  void Sample(const Nucleus& nucleus, const HitNucleonRequest& request, InitialState& out,
              Random& rng) const override;
  bool PauliBlocked(const Nucleus& nucleus, int pdg, const ThreeVector& p) const override;

  // Sample a momentum magnitude from the model's n(k) for Fermi momentum kF.
  double SampleMomentum(double kF, Random& rng) const;
  // Fraction of nucleons above kF and the tail cutoff for Fermi momentum kF.
  double TailFraction(double kF) const;
  double TailCutoff() const { return opt_.bodekRitchie ? kBodekRitchieCutoff : opt_.tailCutoff; }

 private:
  Options opt_;
};

}  // namespace nuflow
