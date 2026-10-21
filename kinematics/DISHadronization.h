#pragma once
// Hadronization of the inclusive DIS hadronic system X (mass W) into a nucleon and pions:
//  - mean charged multiplicity <n_ch> = a + b ln W^2 from deuterium bubble-chamber data
//    (FNAL E545, Zieminska et al., PRD 27, 47 (1983) for nu p / nu n; CERN WA25, Barlag et al.,
//    Z. Phys. C 11, 283 (1982) for nubar p / nubar n);
//  - KNO-type dispersion D_ch = 0.35 <n_ch> (Zieminska et al.);
//  - exact charge and baryon-number conservation, leading baryon from the spectator diquark;
//  - pi0 multiplicity from isospin, <n_pi0> = (n_pi+ + n_pi-)/2;
//  - kinematics: isotropic n-body phase space for W < wLPS, limited-pT longitudinal phase space
//    (rapidity plateau along the momentum transfer, Gaussian pT) above.
// Strange and charmed hadrons are not produced. See theory/dis.md.

#include <utility>
#include <vector>

#include "core/FourVector.h"
#include "core/Random.h"

namespace nuflow {

class DISHadronization {
 public:
  struct Options {
    double dispersionSlope = 0.35;  // D_ch / <n_ch>
    double meanPt2 = 0.25;          // <pT^2> of each hadron w.r.t. the current direction [GeV^2]
    double wLPS = 2.3;              // below: isotropic phase space; above: longitudinal phase space [GeV]
  };
  struct Hadron {
    int pdg;
    FourVector p;
  };

  DISHadronization() = default;
  explicit DISHadronization(Options o) : opt_(o) {}

  // <n_ch>(W) for a CC interaction of `nuPdg` on nucleon `hitPdg` (charged hadrons incl. protons).
  static double MeanChargedMultiplicity(int nuPdg, int hitPdg, double W);

  // Decide the hadron species for a hadronic system of mass W (charge and baryon number of the
  // struck nucleon plus the W boson).
  void SampleSpecies(int nuPdg, int hitPdg, double W, Random& rng, std::vector<int>& pdgs) const;

  // Full hadronization: species + kinematics. X is the hadronic four-momentum (lab frame) and `axis`
  // the direction of the momentum transfer (lab frame). Hadrons are returned in the lab frame.
  void Hadronize(int nuPdg, int hitPdg, const FourVector& X, const ThreeVector& axis, Random& rng,
                 std::vector<Hadron>& out) const;

  // Isotropic n-body phase-space decay of mass W at rest (unweighted). Returns false on failure.
  static bool PhaseSpace(double W, const std::vector<double>& masses, Random& rng, std::vector<FourVector>& p);
  // Longitudinal phase space at rest along +z; the particle with index `backward` gets the
  // smallest rapidity. Returns false if no configuration was found.
  bool LongitudinalPhaseSpace(double W, const std::vector<double>& masses, int backward, Random& rng,
                              std::vector<FourVector>& p) const;

  const Options& Opts() const { return opt_; }

 private:
  Options opt_;
};

}  // namespace nuflow
