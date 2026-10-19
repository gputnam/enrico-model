#pragma once
// NuHepMC v1.0 (HepMC3 ASCII) output, https://github.com/NuHepMC/Spec.
//
// Event graph:
//   [vtx 21: bound-nucleon separation]  target(20) -> struck nucleon(21) + remnant
//   [vtx 1 : primary]                   nu(4) + struck nucleon -> lepton(1) + pre-FSI hadrons(22)
//   [vtx 2 : FSI summary]               pre-FSI hadrons(22) + remnant(23) -> final hadrons(1) + remnant(2009900000, 1)
// For free-nucleon targets the target itself enters the primary vertex.

#include <memory>
#include <string>

#include "orchestration/OutputWriter.h"

namespace HepMC3 {
class GenRunInfo;
class WriterAscii;
}  // namespace HepMC3

namespace nuflow {

class NuHepMCWriter : public OutputWriter {
 public:
  explicit NuHepMCWriter(std::string path);
  ~NuHepMCWriter() override;
  void Begin(const RunInfo& info) override;
  void Write(const Event& event) override;
  void End(const RunSummary& summary) override;

 private:
  std::string path_;
  std::shared_ptr<HepMC3::GenRunInfo> runInfo_;
  std::unique_ptr<HepMC3::WriterAscii> writer_;
};

}  // namespace nuflow
