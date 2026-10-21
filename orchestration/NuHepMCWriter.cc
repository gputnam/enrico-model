#include "orchestration/NuHepMCWriter.h"

#include <HepMC3/Attribute.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenRunInfo.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/WriterAscii.h>

#include <iomanip>
#include <map>
#include <sstream>

#include "core/PDG.h"

namespace nuflow {

namespace {
constexpr int kRemnantPdg = 2009900000;
template <typename A, typename T>
void AddRunAttr(HepMC3::GenRunInfo& r, const std::string& name, const T& value) {
  r.add_attribute(name, std::make_shared<A>(value));
}
HepMC3::FourVector ToHepMC(const FourVector& p) { return {p.px, p.py, p.pz, p.e}; }

// HepMC3's VectorDoubleAttribute prints with std::to_string (6 fixed decimals), which zeroes
// flux contents of O(1e-6) nu/m^2/POT. Write full precision instead; reading is unchanged.
class PreciseVectorDoubleAttribute : public HepMC3::VectorDoubleAttribute {
 public:
  using HepMC3::VectorDoubleAttribute::VectorDoubleAttribute;
  bool to_string(std::string& att) const override {
    std::ostringstream oss;
    oss << std::setprecision(17);
    for (double v : value()) oss << (oss.tellp() > 0 ? " " : "") << v;
    att = oss.str();
    return true;
  }
};
}  // namespace

NuHepMCWriter::NuHepMCWriter(std::string path) : path_(std::move(path)) {}
NuHepMCWriter::~NuHepMCWriter() = default;

void NuHepMCWriter::Begin(const RunInfo& info) {
  using namespace HepMC3;
  runInfo_ = std::make_shared<GenRunInfo>();
  auto& r = *runInfo_;
  r.tools().push_back({"nuflow", "0.1.0", "Staged GeV neutrino event generator"});
  AddRunAttr<IntAttribute>(r, "NuHepMC.Version.Major", 1);
  AddRunAttr<IntAttribute>(r, "NuHepMC.Version.Minor", 0);
  AddRunAttr<IntAttribute>(r, "NuHepMC.Version.Patch", 0);
  AddRunAttr<VectorStringAttribute>(
      r, "NuHepMC.Conventions",
      std::vector<std::string>{"G.C.2", "G.C.4", "E.C.1", "E.C.2", "E.C.3", "V.C.1", "P.C.1", "P.C.2"});
  // Cross sections: 1e-38 cm^2 per target unit; for a single nucleus this is per atom.
  AddRunAttr<DoubleAttribute>(r, "NuHepMC.FluxAveragedTotalCrossSection", info.fluxAveragedSigma);
  AddRunAttr<StringAttribute>(r, "NuHepMC.Units.CrossSection.Unit", "10-38 cm2");
  AddRunAttr<StringAttribute>(r, "NuHepMC.Units.CrossSection.TargetScale",
                              info.target.components.size() == 1 ? "PerAtom" : "PerTargetMolecule");
  r.set_weight_names({"CV"});

  std::vector<int> procIds;
  std::map<int, std::string> procNames;
  for (size_t i = 0; i < info.selector->Size(); ++i) {
    const auto* ch = info.selector->Get(i).channel;
    if (!procNames.count(ch->ProcessId())) procIds.push_back(ch->ProcessId());
    procNames[ch->ProcessId()] = ch->Name();
  }
  AddRunAttr<VectorIntAttribute>(r, "NuHepMC.ProcessIDs", procIds);
  for (const auto& [id, name] : procNames) {
    AddRunAttr<StringAttribute>(r, "NuHepMC.ProcessInfo[" + std::to_string(id) + "].Name", name);
    AddRunAttr<StringAttribute>(r, "NuHepMC.ProcessInfo[" + std::to_string(id) + "].Description",
                                "nuflow " + name);
  }
  AddRunAttr<VectorIntAttribute>(r, "NuHepMC.VertexStatusIDs", std::vector<int>{1, 2, 21});
  AddRunAttr<StringAttribute>(r, "NuHepMC.VertexStatusInfo[1].Name", "PrimaryVertex");
  AddRunAttr<StringAttribute>(r, "NuHepMC.VertexStatusInfo[1].Description", "The primary interaction vertex");
  AddRunAttr<StringAttribute>(r, "NuHepMC.VertexStatusInfo[2].Name", "FSIVertex");
  AddRunAttr<StringAttribute>(r, "NuHepMC.VertexStatusInfo[2].Description", "Final-state interaction summary");
  AddRunAttr<StringAttribute>(r, "NuHepMC.VertexStatusInfo[21].Name", "NucleonSeparationVertex");
  AddRunAttr<StringAttribute>(r, "NuHepMC.VertexStatusInfo[21].Description",
                              "Separation of the struck nucleon from the nuclear remnant");
  AddRunAttr<VectorIntAttribute>(r, "NuHepMC.ParticleStatusIDs", std::vector<int>{1, 4, 20, 21, 22, 23, 25});
  const std::map<int, std::pair<std::string, std::string>> pstat = {
      {1, {"UndecayedPhysical", "Physical final-state particle"}},
      {4, {"IncomingBeam", "Incoming neutrino"}},
      {20, {"Target", "Target nucleus or free nucleon"}},
      {21, {"StruckNucleon", "Bound nucleon struck in the primary interaction"}},
      {22, {"PreFSIHadron", "Hadron leaving the primary vertex, before FSI"}},
      {23, {"NuclearRemnantPreFSI", "Spectator nuclear remnant entering the FSI stage"}},
      {25, {"NuclearRemnantPostFSI", "Excited nuclear remnant after the intranuclear cascade"}}};
  for (const auto& [id, nd] : pstat) {
    AddRunAttr<StringAttribute>(r, "NuHepMC.ParticleStatusInfo[" + std::to_string(id) + "].Name", nd.first);
    AddRunAttr<StringAttribute>(r, "NuHepMC.ParticleStatusInfo[" + std::to_string(id) + "].Description", nd.second);
  }
  std::vector<std::string> abundance;
  for (const auto& c : info.target.components)
    abundance.push_back(std::to_string(c.nucleus.pdg) + "[" + std::to_string(c.count) + "]");
  AddRunAttr<VectorStringAttribute>(r, "NuHepMC.TargetMaterialRelativeAtomicAbundance", abundance);

  // G.C.4 beam description
  const std::string beam = "NuHepMC.Beam[" + std::to_string(info.nuPdg) + "]";
  AddRunAttr<StringAttribute>(r, "NuHepMC.Beam.EnergyUnit", "GEV");
  if (info.flux->IsMonoEnergetic()) {
    AddRunAttr<StringAttribute>(r, beam + ".Type", "MonoEnergetic");
    AddRunAttr<DoubleAttribute>(r, beam + ".MonoEnergetic.Energy", info.flux->EMin());
  } else {
    const auto& h = dynamic_cast<const HistogramFlux&>(*info.flux);
    AddRunAttr<StringAttribute>(r, beam + ".Type", "Histogram");
    AddRunAttr<StringAttribute>(r, "NuHepMC.Beam.FluxUnit", h.Units() + "/GeV");
    AddRunAttr<PreciseVectorDoubleAttribute>(r, beam + ".Histogram.BinEdges", h.Edges());
    AddRunAttr<PreciseVectorDoubleAttribute>(r, beam + ".Histogram.BinContent", h.Density());
    AddRunAttr<BoolAttribute>(r, beam + ".Histogram.ContentIsPerWidth", true);
  }
  writer_ = std::make_unique<WriterAscii>(path_, runInfo_);
}

void NuHepMCWriter::Write(const Event& ev) {
  using namespace HepMC3;
  GenEvent g(runInfo_, Units::GEV, Units::MM);
  g.set_event_number(static_cast<int>(ev.number));
  g.weights().assign(1, ev.weight);
  g.add_attribute("signal_process_id", std::make_shared<IntAttribute>(ev.channel));
  g.add_attribute("lab_pos", std::make_shared<VectorDoubleAttribute>(std::vector<double>{0, 0, 0}));

  const auto& parts = ev.particles;
  std::vector<GenParticlePtr> hp(parts.size());
  auto vPrimary = std::make_shared<GenVertex>();
  vPrimary->set_status(1);
  auto vFSI = std::make_shared<GenVertex>();
  vFSI->set_status(2);
  GenVertexPtr vSep;
  int iTarget = -1, iHit = -1, iRemnant = -1, iRemnantPost = -1;
  for (size_t i = 0; i < parts.size(); ++i) {
    const auto& p = parts[i];
    if (p.status == kStatusRemnantPostFSI) {
      // remnant after the intranuclear cascade: the physical (status-1) remnant of the event
      hp[i] = std::make_shared<GenParticle>(ToHepMC(p.p), kRemnantPdg, 1);
      hp[i]->add_attribute("remnant_particle_number", std::make_shared<IntAttribute>(p.pdg));
      iRemnantPost = static_cast<int>(i);
    } else if (p.status == kStatusRemnant) {
      hp[i] = std::make_shared<GenParticle>(ToHepMC(p.p), kRemnantPdg, 1);
      hp[i]->add_attribute("remnant_particle_number", std::make_shared<IntAttribute>(p.pdg));
      iRemnant = static_cast<int>(i);
    } else {
      hp[i] = std::make_shared<GenParticle>(ToHepMC(p.p), p.pdg, p.status);
    }
    hp[i]->set_generated_mass(p.p.M());
    if (p.status == kStatusTarget) iTarget = static_cast<int>(i);
    if (p.status == kStatusStruckNucleon) iHit = static_cast<int>(i);
  }
  const bool bound = iRemnant >= 0;
  if (bound) {
    vSep = std::make_shared<GenVertex>();
    vSep->set_status(21);
    vSep->add_particle_in(hp[iTarget]);
    vSep->add_particle_out(hp[iHit]);
    vPrimary->add_particle_in(hp[iHit]);
  } else {
    vPrimary->add_particle_in(hp[iTarget]);
  }
  bool anyFSI = false;
  for (size_t i = 0; i < parts.size(); ++i) {
    const auto& p = parts[i];
    switch (p.status) {
      case kStatusBeam: vPrimary->add_particle_in(hp[i]); break;
      case kStatusPreFSI:
        vPrimary->add_particle_out(hp[i]);
        vFSI->add_particle_in(hp[i]);
        anyFSI = true;
        break;
      case kStatusFinal:
        if (p.mother >= 0 && parts[p.mother].status == kStatusPreFSI) vFSI->add_particle_out(hp[i]);
        else vPrimary->add_particle_out(hp[i]);
        break;
      default: break;
    }
  }
  if (bound) {
    // The remnant leaves the separation vertex and passes through the FSI summary vertex.
    // With a cascade the outgoing remnant is the post-FSI one (it balances the FSI vertex).
    auto remnantIn = std::make_shared<GenParticle>(hp[iRemnant]->momentum(), kRemnantPdg, kStatusRemnant);
    remnantIn->set_generated_mass(parts[iRemnant].p.M());
    vSep->add_particle_out(remnantIn);
    vFSI->add_particle_in(remnantIn);
    vFSI->add_particle_out(hp[iRemnantPost >= 0 ? iRemnantPost : iRemnant]);
    g.add_vertex(vSep);
  }
  g.add_vertex(vPrimary);
  if (anyFSI) g.add_vertex(vFSI);

  g.add_attribute("tot_xs", std::make_shared<DoubleAttribute>(ev.totalXSec));
  g.add_attribute("proc_xs", std::make_shared<DoubleAttribute>(ev.processXSec));
  writer_->write_event(g);
}

void NuHepMCWriter::End(const RunSummary&) { writer_->close(); }

}  // namespace nuflow
