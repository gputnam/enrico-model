#pragma once
// Builds the physics components (target, channels, FSI model) from a YAML run configuration.

#include <memory>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "channel/InteractionChannel.h"
#include "core/Nucleus.h"
#include "fsi/FSIModel.h"

namespace nuflow {

// target: [{pdg: 1000060120, count: 1}, {pdg: 1000010010, count: 1}]   or   "1000060120:1,1000010010:1"
Target MakeTarget(const YAML::Node& node);

// channels: [{type: CCQE, nuclear_model: {...}, form_factors: {...}, proposal_q2_scale: 0.6}]
std::vector<std::shared_ptr<InteractionChannel>> MakeChannels(const YAML::Node& node);

// fsi: {type: none}
std::unique_ptr<FSIModel> MakeFSI(const YAML::Node& node);

}  // namespace nuflow
