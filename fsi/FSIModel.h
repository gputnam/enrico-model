#pragma once
// Stage 4: final-state interactions.
//
// An FSIModel takes the hadrons produced at the primary vertex (status kStatusPreFSI) and
// propagates them through the residual nucleus, appending the particles that leave the nucleus
// with status kStatusFinal (and any intermediate particles with generator-specific statuses).

#include <string>

#include "core/Event.h"
#include "core/Nucleus.h"
#include "core/Random.h"

namespace nuflow {

class FSIModel {
 public:
  virtual ~FSIModel() = default;
  virtual std::string Name() const = 0;
  virtual void Propagate(const Nucleus& nucleus, Event& event, Random& rng) const = 0;
};

}  // namespace nuflow
