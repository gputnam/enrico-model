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

#include <memory>
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
  // Pion-multiplicity class of the hadronic final state (theory/dis.md section 3).
  enum class PionClass { kAny, kSingle, kMulti };

  DISHadronization() : DISHadronization(Options{}) {}
  explicit DISHadronization(Options o);

  // <n_ch>(W) for a CC interaction of `nuPdg` on nucleon `hitPdg` (charged hadrons incl. protons).
  static double MeanChargedMultiplicity(int nuPdg, int hitPdg, double W);

  // Decide the hadron species for a hadronic system of mass W (charge and baryon number of the
  // struck nucleon plus the W boson).
  void SampleSpecies(int nuPdg, int hitPdg, double W, Random& rng, std::vector<int>& pdgs) const;
  // As SampleSpecies, restricted to one pion-multiplicity class (exactly one pion, or two or more):
  // the unrestricted multiplicity distribution conditioned on the class (rejection). Returns false if
  // the class is kinematically closed at this W (kMulti below M + 2 m_pi).
  bool SampleSpecies(int nuPdg, int hitPdg, double W, PionClass cls, Random& rng, std::vector<int>& pdgs) const;

  // Probability P_1pi(W) that SampleSpecies produces exactly one pion (nucleon + one pion), tabulated
  // at construction by Monte Carlo with a fixed seed (statistical precision ~0.5%).
  double SinglePionFraction(int nuPdg, int hitPdg, double W) const;

  // Full hadronization: species + kinematics. X is the hadronic four-momentum (lab frame) and `axis`
  // the direction of the momentum transfer (lab frame). Hadrons are returned in the lab frame.
  void Hadronize(int nuPdg, int hitPdg, const FourVector& X, const ThreeVector& axis, Random& rng,
                 std::vector<Hadron>& out, PionClass cls = PionClass::kAny) const;
  // Kinematics only, for a given list of species (first entry: the baryon).
  void HadronizeSpecies(const std::vector<int>& species, const FourVector& X, const ThreeVector& axis, Random& rng,
                        std::vector<Hadron>& out) const;

  // Isotropic n-body phase-space decay of mass W at rest (unweighted). Returns false on failure.
  static bool PhaseSpace(double W, const std::vector<double>& masses, Random& rng, std::vector<FourVector>& p);
  // Longitudinal phase space at rest along +z; the particle with index `backward` gets the
  // smallest rapidity. Returns false if no configuration was found.
  bool LongitudinalPhaseSpace(double W, const std::vector<double>& masses, int backward, Random& rng,
                              std::vector<FourVector>& p) const;

  const Options& Opts() const { return opt_; }

 private:
  struct P1Table {
    double wLo, dW;
    int n;
    std::vector<double> p[4];  // index: 2*(nubar) + (neutron)
  };
  void BuildSinglePionTable();
  void HadronizeImpl(int nuPdg, int hitPdg, const FourVector& X, const ThreeVector& axis, Random& rng,
                     std::vector<Hadron>& out, PionClass cls, const std::vector<int>* fixed) const;
  Options opt_;
  std::shared_ptr<const P1Table> p1_;
};

}  // namespace nuflow
