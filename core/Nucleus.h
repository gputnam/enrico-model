#pragma once
// Nuclear target description and the target mixture.

#include <string>
#include <vector>

namespace nuflow {

struct Nucleus {
  int pdg = 0;
  int Z = 0;
  int A = 0;
  std::string name;
  double mass = 0;          // nuclear (not atomic) mass [GeV]
  double kFermiProton = 0;  // Fermi momenta [GeV]
  double kFermiNeutron = 0;
  double removalEnergy = 0; // mean nucleon removal ("binding") energy for the global Fermi gas [GeV]

  int N() const { return A - Z; }
  bool IsFreeNucleon() const { return A == 1; }
  double KFermi(int nucleonPdg) const;
};

// Parameters per nucleus, read from data/nuclei/nuclei.yaml.
class NuclearDatabase {
 public:
  static const NuclearDatabase& Instance();
  const Nucleus& Get(int pdg) const;
  bool Has(int pdg) const;

 private:
  NuclearDatabase();
  std::vector<Nucleus> nuclei_;
};

// A target is a mixture of nuclei, each with a relative abundance (e.g. CH: C12 x1, H1 x1).
// Cross sections are reported per unit of this mixture (per "molecule").
struct TargetComponent {
  Nucleus nucleus;
  double count = 1;
};

struct Target {
  std::vector<TargetComponent> components;
  std::string name;
  double NucleonsPerUnit() const {
    double n = 0;
    for (const auto& c : components) n += c.count * c.nucleus.A;
    return n;
  }
};

}  // namespace nuflow
