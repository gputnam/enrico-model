#pragma once
// The in-memory event record shared by all four generation stages.

#include <cstdint>
#include <vector>

#include <cmath>

#include "core/FourVector.h"

namespace nuflow {

// Particle status codes follow the NuHepMC conventions where they exist; values >= 21 are
// generator-specific and are declared in the NuHepMC run info.
enum ParticleStatus : int {
  kStatusFinal = 1,           // physical final-state particle (leaves the nucleus)
  kStatusDecayed = 2,
  kStatusBeam = 4,            // incoming neutrino
  kStatusTarget = 20,         // target nucleus
  kStatusStruckNucleon = 21,  // bound nucleon(s) participating in the primary interaction
  kStatusPreFSI = 22,         // hadron produced at the primary vertex, before final-state interactions
  kStatusRemnant = 23,        // spectator nuclear remnant (may be off-shell / excited)
};

struct Particle {
  int pdg = 0;
  int status = 0;
  int mother = -1;  // index into Event::particles; -1 if none
  FourVector p;     // lab-frame four-momentum [GeV]
};

// Initial-state information produced by stage 2 (nuclear model).
struct HitNucleon {
  int pdg = 0;
  ThreeVector p;            // momentum [GeV]
  double energy = 0;        // total (off-shell) energy available to the interaction [GeV]
  double removalEnergy = 0; // energy needed to remove it from the nucleus [GeV]
  double OnShellEnergy(double mass) const { return std::sqrt(p.Mag2() + mass * mass); }
  FourVector Momentum() const { return {energy, p}; }
};

struct InitialState {
  int nHit = 0;
  HitNucleon hit[2];        // QE uses one nucleon, 2p2h will use two
  int remnantPdg = 0;
  FourVector remnant;       // remnant four-momentum = P_A - sum(hit)
};

struct Event {
  int64_t number = 0;
  double weight = 1.0;
  int channel = 0;        // NuHepMC process id of the generated channel
  int nuPdg = 0;
  int targetPdg = 0;
  double enu = 0;
  double totalXSec = 0;    // sigma_tot(E) for this target unit [1e-38 cm^2]
  double processXSec = 0;  // sigma of the selected (channel, nucleus) at E, per nucleus [1e-38 cm^2]
  // Convenience kinematics (true, from the lepton side)
  double Q2 = 0, q0 = 0, q3 = 0, W = 0;
  InitialState initial;
  std::vector<Particle> particles;

  void Clear() {
    particles.clear();
    initial = InitialState{};
    weight = 1.0;
  }
  int Add(int pdg, int status, const FourVector& p, int mother = -1) {
    particles.push_back({pdg, status, mother, p});
    return static_cast<int>(particles.size()) - 1;
  }
};

}  // namespace nuflow
