#pragma once
// Stage 2: initial-state nuclear model.
//
// A NuclearModel samples the bound nucleon(s) that take part in the primary interaction, from the
// model's momentum/removal-energy distribution, independent of the interaction kinematics. The
// interaction-dependent weight (cross section for that initial state) is applied in stage 3, so
// models stay interchangeable between channels.

#include <string>

#include "core/Event.h"
#include "core/FourVector.h"
#include "core/Nucleus.h"
#include "core/Random.h"

namespace nuflow {

// Nuclear medium at the interaction point, for in-medium corrections of the primary vertex (RPA).
struct MediumPoint {
  double rho = 0;       // total nucleon density [fm^-3]
  double rhoCentral = 0;  // density at the nucleus centre [fm^-3] (sets density-dependent couplings)
};

// What an interaction channel needs from the nucleus: how many nucleons, and of which species.
struct HitNucleonRequest {
  int nNucleons = 1;
  int species[2] = {0, 0};  // PDG codes (2212/2112)
};

class NuclearModel {
 public:
  virtual ~NuclearModel() = default;
  virtual std::string Name() const = 0;

  // Fill `out` with the struck nucleon(s) and the spectator remnant.
  virtual void Sample(const Nucleus& nucleus, const HitNucleonRequest& request, InitialState& out,
                      Random& rng) const = 0;

  // True if a nucleon of species `pdg` with momentum `p` cannot be placed in the nucleus.
  virtual bool PauliBlocked(const Nucleus& nucleus, int pdg, const ThreeVector& p) const = 0;
  // The same at the interaction point of `initial` (models with a radial dependence override this).
  virtual bool PauliBlocked(const Nucleus& nucleus, int pdg, const ThreeVector& p, const InitialState& initial) const {
    return PauliBlocked(nucleus, pdg, p);
  }

  // Energy of hit nucleon i entering the primary-vertex two-body kinematics, for an outgoing nucleon of
  // species finalPdg that is produced on shell inside the nucleus. Default: the stage-2 energy, i.e.
  // the outgoing nucleon is asymptotic (global Fermi gas, Smith-Moniz).
  virtual double VertexEnergy(const Nucleus& nucleus, const InitialState& initial, int i, int finalPdg) const {
    return initial.hit[i].energy;
  }
  // Energy an outgoing nucleon of species finalPdg loses when it leaves the nuclear potential: its
  // asymptotic energy is the in-medium (vertex) energy minus this. Default 0.
  virtual double OutgoingPotential(const Nucleus& nucleus, const InitialState& initial, int finalPdg) const { return 0; }

  // Local medium at the interaction point (default: uniform symmetric matter of the global k_F).
  virtual MediumPoint Medium(const Nucleus& nucleus, const InitialState& initial) const;
};

// Helper shared by implementations: fill the remnant after the hit nucleons have been set.
void FillRemnant(const Nucleus& nucleus, InitialState& out);

}  // namespace nuflow
