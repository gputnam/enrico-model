#pragma once
// Flat ROOT output: TTree "events" (one entry per event, particle arrays) plus run metadata:
// TTree "run", sigma(E) TGraphs per process, the flux histogram and the full configuration.

#include <memory>
#include <string>

#include "orchestration/OutputWriter.h"

class TFile;
class TTree;

namespace nuflow {

class RootTreeWriter : public OutputWriter {
 public:
  static constexpr int kMaxParticles = 64;

  // compression: ROOT algorithm*100 + level (e.g. 404 = LZ4 level 4, 505 = ZSTD level 5, 0 = none)
  explicit RootTreeWriter(std::string path, int compression = 505);
  ~RootTreeWriter() override;
  void Begin(const RunInfo& info) override;
  void Write(const Event& event) override;
  void End(const RunSummary& summary) override;

 private:
  std::string path_;
  int compression_;
  std::unique_ptr<TFile> file_;
  TTree* tree_ = nullptr;
  RunInfo info_;

  // Branch buffers
  long long evt_ = 0;
  int channel_ = 0, nuPdg_ = 0, targetPdg_ = 0, hitPdg_ = 0;
  double weight_ = 1, enu_ = 0, Q2_ = 0, q0_ = 0, q3_ = 0, W_ = 0;
  double hitP_[3] = {0, 0, 0};
  double hitE_ = 0, removalE_ = 0;
  int n_ = 0;
  int pdg_[kMaxParticles], status_[kMaxParticles], mother_[kMaxParticles];
  double px_[kMaxParticles], py_[kMaxParticles], pz_[kMaxParticles], E_[kMaxParticles];
};

}  // namespace nuflow
