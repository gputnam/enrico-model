#pragma once
// Stage 4: semi-classical intranuclear cascade (INC) for nucleons and pions
// (theory/fsi_cascade.md).
//
//  * The primary vertex is placed in the nucleus according to the nucleon density
//    (fsi/NuclearDensity); every hadron then moves on a straight line in steps of `step` fm.
//  * Nucleons collide with nucleons of the local Fermi sea (local k_F from the density) with the
//    free pp/np cross sections of Cugnon et al. (elastic, and N N -> N Delta -> N N pi above
//    threshold). Collisions whose outgoing nucleons lie inside the local Fermi sea are Pauli
//    blocked. Nucleons move in a potential of depth V(r) = T_F(r) + S_N (local Fermi kinetic
//    energy plus the measured separation energy); a nucleon whose asymptotic kinetic energy
//    E* = T(r) - V(r) is <= 0 is captured by the remnant.
//  * Pions: Delta(1232)-hole model. Delta formation on a Fermi-sea nucleon with the free P33
//    Breit-Wigner, in-medium width from the Oset-Salcedo Delta self-energy; the Delta either decays
//    (quasi-elastic scattering / charge exchange, Pauli blocked) or is absorbed on one (two-body,
//    pi N N -> N N) or two (three-body) further nucleons. Above the Delta region a non-resonant
//    elastic-like pi N cross section is used.
//  * Hadrons from DIS (process id 600) travel a formation length l = |p| / mu^2 (SKAT) before
//    they can interact.
//  * Escaping hadrons get status 1; the remnant (status 25) takes the remaining four-momentum, so
//    energy, momentum, charge and baryon number are conserved exactly.

#include <vector>

#include "fsi/FSIModel.h"
#include "fsi/NuclearDensity.h"

namespace nuflow {

class IntranuclearCascade : public FSIModel {
 public:
  struct Options {
    double step = 0.2;                     // transport step [fm]
    bool formationZone = true;             // SKAT formation length for DIS hadrons
    double formationMu2 = 0.08;            // mu^2 [GeV^2] (Baranov et al., SKAT)
    std::vector<int> formationProcesses = {600};
    bool nnPionProduction = true;          // N N -> N Delta -> N N pi
    bool pionAbsorption = true;            // Oset-Salcedo two- and three-body absorption
    bool pionQuasiElastic = true;
    bool nucleonRescattering = true;
    double likePairAbsorption = 0.05;      // absorption rate on T=1-only (pp, nn) pairs relative to np pairs
    int maxPauliRetries = 20;              // absorption final states: re-draws before accepting
  };

  // A hadron being transported.
  struct Hadron {
    int pdg = 0;
    ThreeVector x;        // position [fm]
    ThreeVector dir;      // unit direction of motion
    double estar = 0;     // asymptotic kinetic energy [GeV] (nucleons: T_inside - V(r); pions: T)
    double formation = 0; // remaining formation length [fm]
    int ancestor = -1;    // index of the primary (status 22) hadron it descends from
    int generation = 0;   // number of interactions in its history
  };

  // Book-keeping of one cascade (used by the validation tools and tests).
  struct Stats {
    int nucleonCollisions = 0, nucleonPauliBlocked = 0, nnInelastic = 0;
    int pionQE = 0, pionCX = 0, pionAbs2 = 0, pionAbs3 = 0, pionPauliBlocked = 0;
    int captured = 0;
  };

  IntranuclearCascade() = default;
  explicit IntranuclearCascade(const Options& o) : opt_(o) {}
  std::string Name() const override { return "IntranuclearCascade"; }
  void Propagate(const Nucleus& nucleus, Event& event, Random& rng) const override;

  // Transports `active` through the nucleus. On return `escaped` holds the hadrons that left
  // the nucleus (with their exit positions), and zRem / nRem the proton and neutron numbers of
  // the remnant (initially those of the remnant before the cascade).
  void Cascade(const Nucleus& nucleus, const NuclearDensity& density, std::vector<Hadron> active,
               std::vector<Hadron>& escaped, int& zRem, int& nRem, Random& rng, Stats* stats = nullptr) const;

  // Nucleon potential depth V(r) = T_F(r) + S [GeV] for species `pdg`.
  double Potential(const Nucleus& nucleus, const NuclearDensity& density, int pdg, double r) const;

  static const NuclearDensity& DensityFor(const Nucleus& nucleus);  // thread-local cache
  const Options& GetOptions() const { return opt_; }

 private:
  struct Remnant {
    int z, n, z0, n0;
    double Frac(bool proton) const {
      const int have = proton ? z : n, init = proton ? z0 : n0;
      return init > 0 ? std::max(0, have) / double(init) : 0.0;
    }
  };
  // One transport of hadron h; may push secondaries to `active`. Returns when h has escaped,
  // been captured/absorbed, or interacted (then h's products are in `active`).
  void Transport(const Nucleus& nucleus, const NuclearDensity& density, Hadron h, std::vector<Hadron>& active,
                 std::vector<Hadron>& escaped, Remnant& rem, Random& rng, Stats* stats) const;
  bool NucleonInteraction(const Nucleus& nucleus, const NuclearDensity& density, Hadron& h, double r,
                          std::vector<Hadron>& active, Remnant& rem, Random& rng, Stats* stats, bool& consumed) const;
  bool PionInteraction(const Nucleus& nucleus, const NuclearDensity& density, Hadron& h, double r,
                       std::vector<Hadron>& active, Remnant& rem, Random& rng, Stats* stats, bool& consumed) const;
  // Adds a nucleon produced inside the nucleus at position x with in-medium four-momentum p.
  void AddNucleon(const Nucleus& nucleus, const NuclearDensity& density, int pdg, const ThreeVector& x,
                  const FourVector& p, int ancestor, int generation, std::vector<Hadron>& active, Remnant& rem,
                  Stats* stats) const;

  Options opt_;
};

}  // namespace nuflow
