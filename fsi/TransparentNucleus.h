#pragma once
// No final-state interactions: the nucleus is transparent and every primary hadron escapes unchanged.

#include "fsi/FSIModel.h"

namespace nuflow {

class TransparentNucleus : public FSIModel {
 public:
  std::string Name() const override { return "TransparentNucleus"; }
  void Propagate(const Nucleus& nucleus, Event& event, Random& rng) const override;
};

}  // namespace nuflow
