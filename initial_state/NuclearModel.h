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
};

// Helper shared by implementations: fill the remnant after the hit nucleons have been set.
void FillRemnant(const Nucleus& nucleus, InitialState& out);

}  // namespace nuflow
