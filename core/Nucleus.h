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
  double rmsChargeRadius = 0;   // rms charge radius [fm] (sets the density profile used by the cascade)
  double separationProton = 0;  // ground-state proton / neutron separation energies S_p, S_n [GeV]
  double separationNeutron = 0;
  // Ground-state to ground-state Q values of the charged-current transitions [GeV] (nuclear masses):
  // qValueNu = M(Z+1, A) - M(Z, A) (nu: n -> p), qValueNubar = M(Z-1, A) - M(Z, A) (nubar: p -> n).
  // 0 = not given. Used by the local Fermi gas energy balance (theory/local_fermi_gas.md).
  double qValueNu = 0;
  double qValueNubar = 0;

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
