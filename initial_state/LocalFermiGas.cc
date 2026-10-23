#include "initial_state/LocalFermiGas.h"

#include <cmath>
#include <stdexcept>

#include "core/Constants.h"
#include "core/PDG.h"
#include "fsi/NuclearDensity.h"

namespace nuflow {

namespace {

ThreeVector Isotropic(double k, Random& rng) {
  const double c = 2 * rng.Uniform() - 1, s = std::sqrt(std::max(0.0, 1 - c * c));
  const double phi = 2 * constants::kPi * rng.Uniform();
  return {k * s * std::cos(phi), k * s * std::sin(phi), k * c};
}

double FermiEnergy(int pdgCode, double kF) {
  const double m = pdg::NucleonMass(pdgCode);
  return std::sqrt(m * m + kF * kF);
}

}  // namespace

double LocalFermiGas::Well(const Nucleus& nucleus, const NuclearDensity& density, int pdgCode, double r) {
  const double m = pdg::NucleonMass(pdgCode);
  const double kF = density.KFermi(pdgCode, r);
  double S = pdgCode == pdg::kProton ? nucleus.separationProton : nucleus.separationNeutron;
  if (S <= 0) S = 0.008;  // same default as the cascade
  return std::sqrt(kF * kF + m * m) - m + S;
}

double LocalFermiGas::EnergyGap(const Nucleus& nucleus, const NuclearDensity& density, int hitPdg, int finalPdg,
                                double r) {
  if (hitPdg == finalPdg) return 0.0;  // no charge change: ground state -> ground state costs nothing
  const double Q = hitPdg == pdg::kNeutron ? nucleus.qValueNu : nucleus.qValueNubar;
  if (Q == 0) {
    // No tabulated Q value: final nucleon at its Fermi level with zero binding, Delta = S_hit.
    const double S = hitPdg == pdg::kProton ? nucleus.separationProton : nucleus.separationNeutron;
    return S > 0 ? S : 0.008;
  }
  const double eHit = FermiEnergy(hitPdg, density.KFermi(hitPdg, r));
  const double eFin = FermiEnergy(finalPdg, density.KFermi(finalPdg, r));
  return Q - (eFin - eHit);  // Q - Q_LFG(r)
}

void LocalFermiGas::Sample(const Nucleus& nucleus, const HitNucleonRequest& request, InitialState& out,
                           Random& rng) const {
  const NuclearDensity& density = NuclearDensity::For(nucleus);
  out.position = density.SamplePosition(rng);
  out.hasPosition = true;
  const double r = out.position.Mag();
  out.nHit = request.nNucleons;
  for (int i = 0; i < request.nNucleons; ++i) {
    auto& h = out.hit[i];
    h.pdg = request.species[i];
    const double kF = density.KFermi(h.pdg, r);
    h.p = Isotropic(kF * std::cbrt(rng.UniformOpen()), rng);
    const double m = pdg::NucleonMass(h.pdg);
    const double eOn = h.OnShellEnergy(m);
    const double V = Well(nucleus, density, h.pdg, r);
    // Asymptotic-energy convention of the event record: E = sqrt(k^2 + M^2) - V(r), so that the residual
    // nucleus carries the removal energy E_R = V(r) - T(k) as excitation.
    h.removalEnergy = V - (eOn - m);
    h.energy = eOn - V;
  }
  FillRemnant(nucleus, out);
}

bool LocalFermiGas::PauliBlocked(const Nucleus&, int, const ThreeVector&) const {
  throw std::logic_error("LocalFermiGas: Pauli blocking needs the interaction point (use the InitialState overload)");
}

bool LocalFermiGas::PauliBlocked(const Nucleus& nucleus, int pdgCode, const ThreeVector& p,
                                 const InitialState& initial) const {
  const double kF = NuclearDensity::For(nucleus).KFermi(pdgCode, initial.position.Mag());
  return p.Mag2() < kF * kF;
}

double LocalFermiGas::VertexEnergy(const Nucleus& nucleus, const InitialState& initial, int i, int finalPdg) const {
  const auto& h = initial.hit[i];
  if (energy_ == EnergyBalance::kAsymptotic) return h.energy;
  const NuclearDensity& density = NuclearDensity::For(nucleus);
  return h.OnShellEnergy(pdg::NucleonMass(h.pdg)) -
         EnergyGap(nucleus, density, h.pdg, finalPdg, initial.position.Mag());
}

double LocalFermiGas::OutgoingPotential(const Nucleus& nucleus, const InitialState& initial, int finalPdg) const {
  if (energy_ == EnergyBalance::kAsymptotic) return 0.0;
  const NuclearDensity& density = NuclearDensity::For(nucleus);
  const double r = initial.position.Mag();
  const int hitPdg = initial.hit[0].pdg;
  return Well(nucleus, density, hitPdg, r) - EnergyGap(nucleus, density, hitPdg, finalPdg, r);
}

MediumPoint LocalFermiGas::Medium(const Nucleus& nucleus, const InitialState& initial) const {
  const NuclearDensity& density = NuclearDensity::For(nucleus);
  MediumPoint m;
  m.rho = density.Rho(initial.position.Mag());
  m.rhoCentral = density.RhoCentral();
  return m;
}

}  // namespace nuflow
