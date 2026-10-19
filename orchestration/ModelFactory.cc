#include "orchestration/ModelFactory.h"

#include <sstream>
#include <stdexcept>

#include "channel/CCQEChannel.h"
#include "core/Config.h"
#include "fsi/TransparentNucleus.h"
#include "initial_state/FreeNucleon.h"
#include "initial_state/RelativisticFermiGas.h"
#include "kinematics/QuasiElasticKinematics.h"

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
