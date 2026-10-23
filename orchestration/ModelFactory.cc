#include "orchestration/ModelFactory.h"

#include <cmath>
#include <sstream>
#include <stdexcept>

#include "channel/CCQEChannel.h"
#include "channel/CCCOHChannel.h"
#include "channel/CCDISChannel.h"
#include "channel/CCMECChannel.h"
#include "channel/CCRESChannel.h"
#include "channel/CCNStarChannel.h"
#include "core/Config.h"
#include "fsi/IntranuclearCascade.h"
#include "fsi/TransparentNucleus.h"
#include "initial_state/FreeNucleon.h"
#include "initial_state/GroundStateNucleus.h"
#include "initial_state/LocalFermiGas.h"
#include "initial_state/RelativisticFermiGas.h"
#include "initial_state/ShellModelSpectralFunction.h"
#include "kinematics/QuasiElasticKinematics.h"
#include "kinematics/CoherentPionKinematics.h"
#include "kinematics/DISKinematics.h"
#include "kinematics/MECKinematics.h"
#include "kinematics/NStarKinematics.h"
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
  if (type == "SF") {
    ShellModelSpectralFunction::Parameters p;
    p.mfFraction = Get<double>(node, "mf_fraction", p.mfFraction);
    p.e1s = Get<double>(node, "e1s", p.e1s);
    p.e1p = Get<double>(node, "e1p", p.e1p);
    p.widthA = Get<double>(node, "width_a", p.widthA);
    p.widthB = Get<double>(node, "width_b", p.widthB);
    p.eMax = Get<double>(node, "e_max", p.eMax);
    p.kSRC = Get<double>(node, "k_src", p.kSRC);
    p.kMax = Get<double>(node, "k_max", p.kMax);
    p.sigmaCM = Get<double>(node, "sigma_cm", p.sigmaCM);
    p.npFraction = Get<double>(node, "np_fraction", p.npFraction);
    p.emitPartner = Get<bool>(node, "emit_partner", p.emitPartner);
    return std::make_shared<ShellModelSpectralFunction>(p);
  }
  if (type == "LFG") {
    const std::string e = Get<std::string>(node, "energy", "in_medium");
    if (e != "in_medium" && e != "asymptotic") throw std::runtime_error("LFG: energy must be in_medium or asymptotic");
    return std::make_shared<LocalFermiGas>(e == "in_medium" ? LocalFermiGas::EnergyBalance::kInMedium
                                                            : LocalFermiGas::EnergyBalance::kAsymptotic);
  }
  throw std::runtime_error("nuclear_model: unknown type '" + type + "'");
}

// Channels whose kinematics use the global-Fermi-gas energy convention only (the LFG energy balance and
// local Pauli blocking are implemented for CCQE; theory/local_fermi_gas.md).
std::shared_ptr<const NuclearModel> MakeGlobalNuclearModel(const YAML::Node& node, const std::string& channel) {
  const std::string t = Get<std::string>(node, "type", "RFG");
  if (t == "LFG" || t == "SF") throw std::runtime_error(channel + ": nuclear_model " + t + " is only implemented for CCQE");
  return MakeNuclearModel(node);
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
      // Valencia RPA (theory/local_fermi_gas.md): `rpa: true` or a map of parameters.
      const YAML::Node r = c["rpa"];
      if (r && r.IsScalar()) o.rpa = r.as<bool>();
      if (r && r.IsMap()) {
        o.rpa = Get<bool>(r, "enabled", true);
        auto& p = o.rpaParameters;
        p.gPrime = Get<double>(r, "g_prime", p.gPrime);
        p.f2Over4Pi = Get<double>(r, "f2_over_4pi", p.f2Over4Pi);
        p.lambdaPi = Get<double>(r, "lambda_pi", p.lambdaPi);
        p.cRho = Get<double>(r, "c_rho", p.cRho);
        p.lambdaRho = Get<double>(r, "lambda_rho", p.lambdaRho);
        p.mRho = Get<double>(r, "m_rho", p.mRho);
        p.fStarOverF = Get<double>(r, "f_star_over_f", p.fStarOverF);
        p.mDelta = Get<double>(r, "m_delta", p.mDelta);
        p.gammaDelta = Get<double>(r, "gamma_delta", p.gammaDelta);
        p.gammaPH = Get<double>(r, "gamma_ph", p.gammaPH);
        p.c0 = Get<double>(r, "c0", p.c0);
        p.f0pIn = Get<double>(r, "f0p_in", p.f0pIn);
        p.f0pEx = Get<double>(r, "f0p_ex", p.f0pEx);
        p.delta = Get<bool>(r, "delta_hole", p.delta);
        p.qFull = Get<double>(r, "q_full", p.qFull);
        p.qZero = Get<double>(r, "q_zero", p.qZero);
      }
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
      out.push_back(std::make_shared<CCRESChannel>(MakeGlobalNuclearModel(c["nuclear_model"], type), freeNucleon, kin));
    } else if (type == "CCNSTAR") {
      // I = 1/2 resonances of the second resonance region, N(1440) P11, N(1520) D13, N(1535) S11
      // (Lalakulich-Paschos-Piranishvili 2006), W < w_max (= CCRES w_max). theory/nstar.md.
      NStarKinematics::Options o;
      o.wMax = Get<double>(c, "w_max", o.wMax);
      o.MV = Get<double>(c, "MV", o.MV);
      o.MA = Get<double>(c, "MA", o.MA);
      o.proposalQ2Scale = Get<double>(c, "proposal_q2_scale", o.proposalQ2Scale);
      if (c["resonances"]) {
        o.resonances.clear();
        for (const auto& r : c["resonances"]) o.resonances.push_back(r.as<std::string>());
      }
      auto kin = std::make_shared<NStarKinematics>(o);
      out.push_back(std::make_shared<CCNStarChannel>(MakeGlobalNuclearModel(c["nuclear_model"], type), freeNucleon, kin));
    } else if (type == "CCDIS") {
      // Inelastic continuum above w_min from the Bodek-Yang structure functions (theory/dis.md).
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
      o.wMax = Get<double>(c, "w_max", 0.0);
      const YAML::Node sp = c["single_pion"];
      o.scaleI32 = Get<double>(sp, "scale_i32", 1.0);
      o.scaleI12 = Get<double>(sp, "scale_i12", 1.0);
      o.isospinDecomposition = Get<bool>(sp, "isospin", false);
      o.proposalQ2Scale = Get<double>(c, "proposal_q2_scale", o.proposalQ2Scale);
      if (o.scaleI32 < 0 || o.scaleI12 < 0) throw std::runtime_error(type + ": single_pion scales must be >= 0");
      auto kin = std::make_shared<DISKinematics>(BodekYang(bp), DISHadronization(ho), o);
      out.push_back(std::make_shared<CCDISChannel>(MakeGlobalNuclearModel(c["nuclear_model"], type), freeNucleon, kin));
    } else if (type == "CCMEC") {
      // 2p2h meson-exchange currents, relativistic MEC of Ruiz Simo et al. (theory/mec.md).
      MECParameters par;
      const YAML::Node m = c["mec"];
      par.C5A0 = Get<double>(m, "C5A0", par.C5A0);
      par.MADelta = Get<double>(m, "MA", par.MADelta);
      par.C3V0 = Get<double>(m, "C3V0", par.C3V0);
      par.MV = Get<double>(m, "MV", par.MV);
      par.fStar = Get<double>(m, "f_star", par.fStar);
      par.MDelta = Get<double>(m, "mass", par.MDelta);
      par.gA = Get<double>(m, "gA", par.gA);
      par.fPi = Get<double>(m, "f_pi", par.fPi);
      par.realDeltaPropagator = Get<std::string>(m, "delta_propagator", "real") == "real";
      par.lambdaPiNN = Get<double>(m, "lambda_piNN", par.lambdaPiNN);
      par.lambdaPiNDelta = Get<double>(m, "lambda_piNDelta", par.lambdaPiNDelta);
      // Nucleon form factors (only the isovector Dirac F1V enters): HNV dipole/Galster by default.
      YAML::Node nff = c["nucleon_form_factors"] ? YAML::Clone(c["nucleon_form_factors"]) : YAML::Node(YAML::NodeType::Map);
      if (!nff["vector"]) nff["vector"] = "dipole";
      if (!nff["axial"]) nff["axial"] = "dipole";
      if (!nff["gA"]) nff["gA"] = par.gA;
      if (!nff["MA"]) nff["MA"] = 1.05;
      MECKinematics::Options o;
      o.proposalQScale = Get<double>(c, "proposal_q_scale", o.proposalQScale);
      o.q3Max = Get<double>(c, "q3_max", o.q3Max);
      auto kin = std::make_shared<MECKinematics>(MECCurrent(par, NucleonFormFactors::FromConfig(nff)), o);
      out.push_back(std::make_shared<CCMECChannel>(MakeGlobalNuclearModel(c["nuclear_model"], type), freeNucleon, kin));
    } else if (type == "CCCOH") {
      // CC coherent pion production, Berger-Sehgal PCAC model with pi-12C elastic data and the
      // lepton-mass (pion-pole) correction (theory/coherent_pion.md). No nuclear_model: the nucleus
      // stays in its ground state.
      PionNucleusElastic::Parameters ep;
      const YAML::Node e = c["elastic"];
      ep.R0 = Get<double>(e, "R0", ep.R0);
      ep.tableMaxT = Get<double>(e, "table_max_T", ep.tableMaxT);
      ep.rsMinT = Get<double>(e, "rs_min_T", ep.rsMinT);
      if (ep.rsMinT < ep.tableMaxT) throw std::runtime_error(type + ": rs_min_T must be >= table_max_T");
      CoherentPionKinematics::Options o;
      o.fPi = Get<double>(c, "f_pi", o.fPi);
      o.mA = Get<double>(c, "mA", o.mA);
      o.leptonMassCorrection = Get<bool>(c, "lepton_mass_correction", o.leptonMassCorrection);
      o.proposalSlope = Get<double>(c, "proposal_slope", o.proposalSlope);
      auto kin = std::make_shared<CoherentPionKinematics>(PionNucleusElastic(ep), o);
      out.push_back(std::make_shared<CCCOHChannel>(std::make_shared<GroundStateNucleus>(), freeNucleon, kin));
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
  if (type == "cascade") {
    // Semi-classical intranuclear cascade for nucleons and pions (theory/fsi_cascade.md).
    IntranuclearCascade::Options o;
    o.step = Get<double>(node, "step", o.step);
    o.nucleonRescattering = Get<bool>(node, "nucleon_rescattering", o.nucleonRescattering);
    o.nnPionProduction = Get<bool>(node, "nn_pion_production", o.nnPionProduction);
    o.pionQuasiElastic = Get<bool>(node, "pion_quasielastic", o.pionQuasiElastic);
    o.pionAbsorption = Get<bool>(node, "pion_absorption", o.pionAbsorption);
    o.likePairAbsorption = Get<double>(node, "like_pair_absorption", o.likePairAbsorption);
    o.formationZone = Get<bool>(node, "formation_zone", o.formationZone);
    o.formationMu2 = Get<double>(node, "formation_mu2", o.formationMu2);
    if (node && node["formation_processes"]) o.formationProcesses = node["formation_processes"].as<std::vector<int>>();
    return std::make_unique<IntranuclearCascade>(o);
  }
  throw std::runtime_error("fsi: unknown type '" + type + "'");
}

}  // namespace nuflow
