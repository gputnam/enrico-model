#include "orchestration/RootTreeWriter.h"

#include <algorithm>
#include <stdexcept>

#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TTree.h>

namespace nuflow {

RootTreeWriter::RootTreeWriter(std::string path, int compression)
    : path_(std::move(path)), compression_(compression) {}
RootTreeWriter::~RootTreeWriter() = default;

void RootTreeWriter::Begin(const RunInfo& info) {
  info_ = info;
  // Let ROOT compress baskets in parallel (TTree::Fill flushes clusters through its task pool);
  // otherwise the single writer thread is the bottleneck of the whole generator.
  if (!ROOT::IsImplicitMTEnabled()) ROOT::EnableImplicitMT(4);
  file_.reset(TFile::Open(path_.c_str(), "RECREATE"));
  if (!file_ || file_->IsZombie()) throw std::runtime_error("RootTreeWriter: cannot create " + path_);
  file_->SetCompressionSettings(compression_);
  tree_ = new TTree("events", "nuflow events");
  tree_->SetAutoFlush(-32000000);  // ~32 MB clusters: fewer, larger baskets
  tree_->Branch("evt", &evt_, "evt/L");
  tree_->Branch("weight", &weight_, "weight/D");
  tree_->Branch("channel", &channel_, "channel/I");
  tree_->Branch("nu_pdg", &nuPdg_, "nu_pdg/I");
  tree_->Branch("target_pdg", &targetPdg_, "target_pdg/I");
  tree_->Branch("Enu", &enu_, "Enu/D");
  tree_->Branch("Q2", &Q2_, "Q2/D");
  tree_->Branch("q0", &q0_, "q0/D");
  tree_->Branch("q3", &q3_, "q3/D");
  tree_->Branch("W", &W_, "W/D");
  tree_->Branch("hit_pdg", &hitPdg_, "hit_pdg/I");
  tree_->Branch("hit_p", hitP_, "hit_p[3]/D");
  tree_->Branch("hit_E", &hitE_, "hit_E/D");
  tree_->Branch("removal_E", &removalE_, "removal_E/D");
  tree_->Branch("n", &n_, "n/I");
  tree_->Branch("pdg", pdg_, "pdg[n]/I");
  tree_->Branch("status", status_, "status[n]/I");
  tree_->Branch("mother", mother_, "mother[n]/I");
  tree_->Branch("E", E_, "E[n]/D");
  tree_->Branch("px", px_, "px[n]/D");
  tree_->Branch("py", py_, "py[n]/D");
  tree_->Branch("pz", pz_, "pz[n]/D");
}

void RootTreeWriter::Write(const Event& ev) {
  evt_ = ev.number;
  weight_ = ev.weight;
  channel_ = ev.channel;
  nuPdg_ = ev.nuPdg;
  targetPdg_ = ev.targetPdg;
  enu_ = ev.enu;
  Q2_ = ev.Q2;
  q0_ = ev.q0;
  q3_ = ev.q3;
  W_ = ev.W;
  const auto& h = ev.initial.hit[0];
  hitPdg_ = h.pdg;
  hitP_[0] = h.p.x;
  hitP_[1] = h.p.y;
  hitP_[2] = h.p.z;
  hitE_ = h.energy;
  removalE_ = h.removalEnergy;
  n_ = static_cast<int>(std::min<size_t>(ev.particles.size(), kMaxParticles));
  for (int i = 0; i < n_; ++i) {
    const auto& p = ev.particles[i];
    pdg_[i] = p.pdg;
    status_[i] = p.status;
    mother_[i] = p.mother;
    E_[i] = p.p.e;
    px_[i] = p.p.px;
    py_[i] = p.p.py;
    pz_[i] = p.p.pz;
  }
  tree_->Fill();
}

void RootTreeWriter::End(const RunSummary& s) {
  file_->cd();
  tree_->Write();

  TTree run("run", "nuflow run summary");
  long long nEvents = s.nEvents, nTrials = s.nTrials, nOver = s.nOverweight;
  double fluxSigma = info_.fluxAveragedSigma, fluxInt = info_.fluxIntegral, maxOver = s.maxOverweightRatio;
  double nucleons = info_.target.NucleonsPerUnit(), wall = s.wallSeconds;
  int nu = info_.nuPdg;
  unsigned long long seed = info_.seed;
  run.Branch("n_events", &nEvents, "n_events/L");
  run.Branch("n_trials", &nTrials, "n_trials/L");
  run.Branch("n_overweight", &nOver, "n_overweight/L");
  run.Branch("max_overweight_ratio", &maxOver, "max_overweight_ratio/D");
  run.Branch("flux_avg_xsec", &fluxSigma, "flux_avg_xsec/D");  // 1e-38 cm^2 per target unit
  run.Branch("nucleons_per_target", &nucleons, "nucleons_per_target/D");
  run.Branch("flux_integral", &fluxInt, "flux_integral/D");
  run.Branch("nu_pdg", &nu, "nu_pdg/I");
  run.Branch("seed", &seed, "seed/l");
  run.Branch("wall_seconds", &wall, "wall_seconds/D");
  run.Fill();
  run.Write();

  TNamed("config", info_.configYaml.c_str()).Write();
  TNamed("target", info_.target.name.c_str()).Write();
  TNamed("flux_description", info_.flux->Description().c_str()).Write();
  TNamed("fsi_model", info_.fsiModel.c_str()).Write();

  for (size_t i = 0; i < info_.selector->Size(); ++i) {
    const auto& p = info_.selector->Get(i);
    const auto& t = p.table;
    TGraphErrors g(static_cast<int>(t.energy.size()));
    for (size_t j = 0; j < t.energy.size(); ++j) {
      g.SetPoint(static_cast<int>(j), t.energy[j], t.sigma[j]);
      g.SetPointError(static_cast<int>(j), 0, t.sigmaErr[j]);
    }
    const std::string name = "xsec_" + p.channel->Name() + "_" + p.nucleus.name;
    g.SetName(name.c_str());
    g.SetTitle((p.channel->Name() + " on " + p.nucleus.name + ";E_{#nu} [GeV];#sigma [10^{-38} cm^{2}/nucleus]").c_str());
    g.Write();
  }
  if (const auto* h = dynamic_cast<const HistogramFlux*>(info_.flux)) {
    const auto& e = h->Edges();
    TH1D flux("flux", (h->Description() + ";E_{#nu} [GeV];flux density [" + h->Units() + "/GeV]").c_str(),
              static_cast<int>(h->NBins()), e.data());
    for (size_t i = 0; i < h->NBins(); ++i) flux.SetBinContent(static_cast<int>(i) + 1, h->Density()[i]);
    flux.Write();
  }
  file_->Close();
}

}  // namespace nuflow
