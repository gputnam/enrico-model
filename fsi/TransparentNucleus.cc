#include "fsi/TransparentNucleus.h"

namespace nuflow {

void TransparentNucleus::Propagate(const Nucleus&, Event& event, Random&) const {
  const size_t n = event.particles.size();
  for (size_t i = 0; i < n; ++i) {
    const Particle part = event.particles[i];  // copy: Add() may reallocate
    if (part.status == kStatusPreFSI) event.Add(part.pdg, kStatusFinal, part.p, static_cast<int>(i));
  }
}

}  // namespace nuflow
