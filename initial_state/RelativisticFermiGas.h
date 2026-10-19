#pragma once
// Global relativistic Fermi gas (Smith & Moniz, Nucl. Phys. B43 (1972) 605), with an optional
// Bodek-Ritchie-style high-momentum tail n(k) ~ 1/k^4 above k_F (Phys. Rev. D 23, 1070 (1981)).

#include "initial_state/NuclearModel.h"

namespace nuflow {

class RelativisticFermiGas : public NuclearModel {
 public:
  struct Options {
    double tailFraction = 0.0;  // fraction of nucleons in the k > k_F tail (0 = pure Smith-Moniz)
    double tailCutoff = 0.5;    // maximum tail momentum [GeV]
  };

  RelativisticFermiGas() = default;
  explicit RelativisticFermiGas(Options opt) : opt_(opt) {}

  std::string Name() const override { return opt_.tailFraction > 0 ? "RFG+tail" : "RFG"; }
  void Sample(const Nucleus& nucleus, const HitNucleonRequest& request, InitialState& out,
              Random& rng) const override;
  bool PauliBlocked(const Nucleus& nucleus, int pdg, const ThreeVector& p) const override;

  // Sample a momentum magnitude from the model's n(k) for Fermi momentum kF.
  double SampleMomentum(double kF, Random& rng) const;

 private:
  Options opt_;
};

}  // namespace nuflow
