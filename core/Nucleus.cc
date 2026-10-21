#include "core/Nucleus.h"

#include <stdexcept>

#include "core/Config.h"
#include "core/PDG.h"

namespace nuflow {

double Nucleus::KFermi(int nucleonPdg) const {
  return nucleonPdg == pdg::kProton ? kFermiProton : kFermiNeutron;
}

const NuclearDatabase& NuclearDatabase::Instance() {
  static const NuclearDatabase db;
  return db;
}

NuclearDatabase::NuclearDatabase() {
  const YAML::Node root = LoadConfig("data/nuclei/nuclei.yaml");
  for (const auto& n : root["nuclei"]) {
    Nucleus nuc;
    nuc.name = n["name"].as<std::string>();
    nuc.Z = n["Z"].as<int>();
    nuc.A = n["A"].as<int>();
    nuc.pdg = pdg::IonCode(nuc.Z, nuc.A);
    nuc.mass = n["mass"].as<double>();
    const double kF = nuflow::Get<double>(n, "kF", 0.0);
    nuc.kFermiProton = nuflow::Get<double>(n, "kF_p", kF);
    nuc.kFermiNeutron = nuflow::Get<double>(n, "kF_n", kF);
    nuc.removalEnergy = nuflow::Get<double>(n, "Eb", 0.0);
    nuc.rmsChargeRadius = nuflow::Get<double>(n, "rch", 0.0);
    nuc.separationProton = nuflow::Get<double>(n, "Sp", 0.0);
    nuc.separationNeutron = nuflow::Get<double>(n, "Sn", 0.0);
    nuclei_.push_back(nuc);
  }
}

// Accept both 2212 and 1000010010 for hydrogen (and 2112 / 1000000010 for a free neutron).
static int Canonical(int code) { return pdg::IonCode(pdg::IonZ(code), pdg::IonA(code)); }

bool NuclearDatabase::Has(int code) const {
  code = Canonical(code);
  for (const auto& n : nuclei_)
    if (n.pdg == code) return true;
  return false;
}

const Nucleus& NuclearDatabase::Get(int code) const {
  code = Canonical(code);
  for (const auto& n : nuclei_)
    if (n.pdg == code) return n;
  throw std::runtime_error("nuflow: nucleus " + std::to_string(code) + " not in data/nuclei/nuclei.yaml");
}

}  // namespace nuflow
