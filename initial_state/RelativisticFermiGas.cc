#include "initial_state/RelativisticFermiGas.h"

#include <cmath>

#include "core/Constants.h"
#include "core/PDG.h"

namespace nuflow {

double RelativisticFermiGas::SampleMomentum(double kF, Random& rng) const {
  // Fermi sea: n(k) ~ k^2 on [0, kF]  ->  k = kF * u^(1/3).
  if (opt_.tailFraction <= 0 || rng.Uniform() >= opt_.tailFraction) return kF * std::cbrt(rng.UniformOpen());
  // Tail: k^2 n(k) ~ k^2 / k^4 = 1/k^2 on [kF, kc]  ->  inverse CDF of 1/k^2.
  const double kc = opt_.tailCutoff;
  const double u = rng.Uniform();
  return 1.0 / (1.0 / kF - u * (1.0 / kF - 1.0 / kc));
}

void RelativisticFermiGas::Sample(const Nucleus& nucleus, const HitNucleonRequest& request, InitialState& out,
                                  Random& rng) const {
  out.nHit = request.nNucleons;
  for (int i = 0; i < request.nNucleons; ++i) {
    auto& h = out.hit[i];
    h.pdg = request.species[i];
    const double k = SampleMomentum(nucleus.KFermi(h.pdg), rng);
    const double cosT = 2 * rng.Uniform() - 1;
    const double sinT = std::sqrt(std::max(0.0, 1 - cosT * cosT));
    const double phi = 2 * constants::kPi * rng.Uniform();
    h.p = {k * sinT * std::cos(phi), k * sinT * std::sin(phi), k * cosT};
    // Smith-Moniz: the bound nucleon has energy sqrt(k^2 + M^2) - epsilon_b.
    h.removalEnergy = nucleus.removalEnergy;
    h.energy = h.OnShellEnergy(pdg::NucleonMass(h.pdg)) - h.removalEnergy;
  }
  FillRemnant(nucleus, out);
}

bool RelativisticFermiGas::PauliBlocked(const Nucleus& nucleus, int pdg, const ThreeVector& p) const {
  const double kF = nucleus.KFermi(pdg);
  return p.Mag2() < kF * kF;
}

}  // namespace nuflow
