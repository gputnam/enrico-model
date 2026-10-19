#pragma once
// Event output interface. Writers receive events in a deterministic order from a single thread.

#include <string>
#include <vector>

#include "channel/ChannelSelector.h"
#include "core/Event.h"
#include "core/Nucleus.h"
#include "flux/FluxSource.h"

namespace nuflow {

struct RunInfo {
  std::string configYaml;
  int nuPdg = 0;
  Target target;
  const FluxSource* flux = nullptr;
  const ChannelSelector* selector = nullptr;
  std::string fsiModel;
  uint64_t seed = 0;
  double fluxAveragedSigma = 0;  // per target unit [1e-38 cm^2]
  double fluxIntegral = 0;       // flux file units (0 for mono-energetic)
  std::string fluxUnits;
};

struct RunSummary {
  long long nEvents = 0;
  long long nTrials = 0;
  long long nOverweight = 0;
  double maxOverweightRatio = 0;  // max w / w_envelope over all trials (> 1 means envelope violated)
  double wallSeconds = 0;
};

class OutputWriter {
 public:
  virtual ~OutputWriter() = default;
  virtual void Begin(const RunInfo& info) = 0;
  virtual void Write(const Event& event) = 0;
  virtual void End(const RunSummary& summary) = 0;
};

}  // namespace nuflow
