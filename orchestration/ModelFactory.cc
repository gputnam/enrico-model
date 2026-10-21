#include "orchestration/ModelFactory.h"

#include <sstream>
#include <stdexcept>

#include "channel/CCQEChannel.h"
#include "channel/CCDISChannel.h"
#include "channel/CCRESChannel.h"
#include "core/Config.h"
#include "fsi/TransparentNucleus.h"
#include "initial_state/FreeNucleon.h"
#include "initial_state/RelativisticFermiGas.h"
#include "kinematics/QuasiElasticKinematics.h"
#include "kinematics/DISKinematics.h"
#include "kinematics/ResonantPionKinematics.h"

namespace nuflow {

Target MakeTarget(const YAML::Node& node) {
  Target t;
  const auto& db = NuclearDatabase::Instance();
  auto add = [&](int code, double count) {
    t.components.push_back({db.Get(code), count});
    if (!t.name.empty()) t.name += " + ";
    std::ostringstream s;
    s << count << "x" << db.Get(code).name;
    t.name += s.str();
  };
  if (node.IsScalar()) {
    std::string spec = node.as<std::string>();
    std::stringstream ss(spec);
    for (std::string item; std::getline(ss, item, ',');) {
      const auto colon = item.find(':');
      add(std::stoi(item.substr(0, colon)), colon == std::string::npos ? 1.0 : std::stod(item.substr(colon + 1)));
    }
  } else {
    for (const auto& c : node) add(c["pdg"].as<int>(), Get<double>(c, "count", 1.0));
  }
  if (t.components.empty()) throw std::runtime_error("target: no components");
  return t;
}

namespace {

std::shared_ptr<const NuclearModel> MakeNuclearModel(const YAML::Node& node) {
  const std::string type = Get<std::string>(node, "type", "RFG");
  if (type == "RFG") {
    RelativisticFermiGas::Options o;
    o.bodekRitchie = Get<std::string>(node, "tail", "none") == "bodek_ritchie";
    o.tailFraction = Get<double>(node, "tail_fraction", 0.0);
    o.tailCutoff = Get<double>(node, "tail_cutoff", 0.5);
    return std::make_shared<RelativisticFermiGas>(o);
  }
  throw std::runtime_error("nuclear_model: unknown type '" + type + "'");
}

}  // namespace

std::vector<std::shared_ptr<InteractionChannel>> MakeChannels(const YAML::Node& node) {
  std::vector<std::shared_ptr<InteractionChannel>> out;
  auto freeNucleon = std::make_shared<FreeNucleon>();
  for (const auto& c : node) {
    const std::string type = c["type"].as<std::string>();
    if (type == "CCQE") {
      QuasiElasticKinematics::Options o;
      o.proposalQ2Scale = Get<double>(c, "proposal_q2_scale", o.proposalQ2Scale);
      auto kin = std::make_shared<QuasiElasticKinematics>(NucleonFormFactors::FromConfig(c["form_factors"]), o);
      out.push_back(std::make_shared<CCQEChannel>(MakeNuclearModel(c["nuclear_model"]), freeNucleon, kin));
    } else if (type == "CCRES") {
      // Delta(1232) region pion production, HNV model (theory/pion_production.md).
      HNVParameters par;
      const YAML::Node d = c["delta"];
      par.C5A0 = Get<double>(d, "C5A0", par.C5A0);
      par.MADelta = Get<double>(d, "MA", par.MADelta);
      par.MDelta = Get<double>(d, "mass", par.MDelta);
      par.fStar = Get<double>(d, "f_star", par.fStar);
      par.background = Get<bool>(c, "background", par.background);
      par.gA = Get<double>(c["nucleon_form_factors"], "gA", par.gA);
      // Nucleon form factors of the background terms (HNV: dipole/Galster vector, g_A = 1.26, M_A = 1.05).
      YAML::Node nff = c["nucleon_form_factors"] ? YAML::Clone(c["nucleon_form_factors"]) : YAML::Node(YAML::NodeType::Map);
      if (!nff["vector"]) nff["vector"] = "dipole";
      if (!nff["axial"]) nff["axial"] = "dipole";
      if (!nff["gA"]) nff["gA"] = par.gA;
      if (!nff["MA"]) nff["MA"] = 1.05;
      ResonantPionKinematics::Options o;
      o.wMax = Get<double>(c, "w_max", o.wMax);
      o.proposalQ2Scale = Get<double>(c, "proposal_q2_scale", o.proposalQ2Scale);
      o.bwFraction = Get<double>(c, "proposal_bw_fraction", o.bwFraction);
      auto kin = std::make_shared<ResonantPionKinematics>(
          PionProductionCurrent(par, NucleonFormFactors::FromConfig(nff)), o);
      out.push_back(std::make_shared<CCRESChannel>(MakeNuclearModel(c["nuclear_model"]), freeNucleon, kin));
    } else if (type == "CCDIS") {
      // Inclusive inelastic scattering above w_min: Bodek-Yang structure functions (theory/dis.md).
      BodekYang::Parameters bp;
      const YAML::Node b = c["bodek_yang"];
      bp.A = Get<double>(b, "A", bp.A);
      bp.B = Get<double>(b, "B", bp.B);
      bp.N = Get<double>(b, "N", bp.N);
      bp.axialTypeII = Get<std::string>(b, "axial", "typeII") == "typeII";
      bp.charm = Get<bool>(b, "charm", bp.charm);
      bp.duCorrection = Get<bool>(b, "du_correction", bp.duCorrection);
      DISHadronization::Options ho;
      const YAML::Node h = c["hadronization"];
      ho.dispersionSlope = Get<double>(h, "dispersion_slope", ho.dispersionSlope);
      ho.meanPt2 = Get<double>(h, "mean_pt2", ho.meanPt2);
      ho.wLPS = Get<double>(h, "w_lps", ho.wLPS);
      DISKinematics::Options o;
      o.wMin = Get<double>(c, "w_min", o.wMin);
      o.proposalQ2Scale = Get<double>(c, "proposal_q2_scale", o.proposalQ2Scale);
      auto kin = std::make_shared<DISKinematics>(BodekYang(bp), DISHadronization(ho), o);
      out.push_back(std::make_shared<CCDISChannel>(MakeNuclearModel(c["nuclear_model"]), freeNucleon, kin));
    } else {
      throw std::runtime_error("channels: unknown type '" + type + "'");
    }
  }
  if (out.empty()) throw std::runtime_error("channels: none configured");
  return out;
}

std::unique_ptr<FSIModel> MakeFSI(const YAML::Node& node) {
  const std::string type = Get<std::string>(node, "type", "none");
  if (type == "none") return std::make_unique<TransparentNucleus>();
  throw std::runtime_error("fsi: unknown type '" + type + "'");
}

}  // namespace nuflow
