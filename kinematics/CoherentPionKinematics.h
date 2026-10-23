#pragma once
// Charged-current coherent pion production, nu A -> l- pi+ A (nubar A -> l+ pi- A), in the PCAC model of
// Ch. Berger and L.M. Sehgal, Phys. Rev. D 79, 053003 (2009) [arXiv:0812.2653] (BS09):
//
//   dsigma/dQ2 dy dt = G_F^2 cos^2(theta_C) f_pi^2 / (2 pi^2) * (E/|q|) u v
//                      * [ (G_A - Q2min/(2(Q2+m_pi^2)))^2 + y/4 (Q2 - Q2min) Q2min/(Q2+m_pi^2)^2 ]
//                      * dsigma(pi A -> pi A)/dt (T_pi, |t|)                                 (BS09 Eq. 6)
//
// with u, v = (E + E' +- |q|)/2E, G_A = m_A^2/(m_A^2 + Q2), Q2min = m_l^2 y/(1-y) (Adler / Kopeliovich-
// Marage lepton-mass correction with the pion pole). The pion-nucleus elastic cross section is
// dsigma/dt = A1(T_pi) exp(-b1(T_pi)|t|): for 12C below T_pi = 0.87 GeV, the BS09 fit to the measured
// pi+- 12C elastic differential cross sections (phase-shift analyses of the pion-carbon data, BS09 Fig. 2
// table); above, a Glauber multiple-scattering amplitude (GlauberPionNucleus; default since v8) or, with
// elastic.model = rein_sehgal, the Rein-Sehgal hadronic model (NPB 223, 29 (1983); BS09 Eqs. 10-15) above
// 1.7 GeV with a linear (A1, b1) bridge (v7). The final state is
// l + pi + A in its ground state with exact two-body recoil kinematics. See theory/coherent_pion.md.

#include <string>

#include "kinematics/GlauberPionNucleus.h"
#include "kinematics/KinematicsGenerator.h"

namespace nuflow {

// Pion-nucleus elastic differential cross section (isoscalar, pi+ = pi-).
//  - T <= tableMaxT: BS09 12C data fit A1 exp(-b1 |t|) (other A: A^(2/3) scaling).
//  - model "glauber" (default): Glauber amplitude (GlauberPionNucleus) above tableMaxT, blended linearly in
//    dsigma/dt with the last data node over [tableMaxT, tableMaxT + blendWidth].
//  - model "rein_sehgal" (v7): linear (A1, b1) bridge to the Rein-Sehgal optical model at rsMinT.
class PionNucleusElastic {
 public:
  struct Parameters {
    std::string model = "glauber";  // "glauber" | "rein_sehgal"
    double R0 = 1.057;          // nuclear radius parameter of the RS optical model [fm] (BS09 Eq. 12)
    double tableMaxT = 0.870;   // last T_pi of the 12C data fit [GeV]
    double rsMinT = 1.7;        // T_pi above which the RS hadronic model is used [GeV] (BS09 section 5)
    double blendWidth = 0.1;    // glauber: data -> Glauber blending interval [GeV]
  };
  explicit PionNucleusElastic(Parameters p) : p_(p) {}
  PionNucleusElastic() : PionNucleusElastic(Parameters{}) {}

  // A1 [mb/GeV^2] and b1 [GeV^-2] at pion kinetic energy T [GeV] for a nucleus of mass number A.
  void Coefficients(double T, int A, double& A1, double& b1) const;
  // dsigma/dt [mb/GeV^2].
  double DSigmaDt(double T, double absT, int A) const;
  // Integrated elastic cross section [mb] (validation; A1/b1 where exponential).
  double Sigma(double T, int A) const;
  bool UsesGlauber() const { return p_.model == "glauber"; }
  const GlauberPionNucleus& Glauber() const { return glauber_; }

  // Rein-Sehgal hadronic model pieces (BS09 Eqs. 10-15), public for validation.
  static double PiNTotal(double T);    // (sigma(pi+ p) + sigma(pi- p))/2 [mb], PDG Regge fit
  static double PiNElastic(double T);  // isoscalar pi N elastic [mb], diffraction-peak estimate
  void ReinSehgal(double T, int A, double& A1, double& b1) const;

  const Parameters& Params() const { return p_; }

 private:
  Parameters p_;
  GlauberPionNucleus glauber_;
};

class CoherentPionKinematics : public KinematicsGenerator {
 public:
  struct Options {
    double fPi = 0.1307;        // pion decay constant [GeV] (BS09: 130.7 MeV)
    double mA = 0.95;           // axial-vector form factor mass [GeV] (BS09 Eq. 5)
    bool leptonMassCorrection = true;  // pion-pole (Adler) term of BS09 Eq. 6
    double proposalSlope = 0.8;  // |t| proposal slope / b_RS(A)
  };

  CoherentPionKinematics(PionNucleusElastic elastic, Options opt) : elastic_(elastic), opt_(opt) {}

  std::string Name() const override;
  double Generate(int nuPdg, double enu, const Nucleus& nucleus, const InitialState& initial,
                  const NuclearModel& nuclearModel, Event& event, Random& rng) const override;

  // dsigma/(dQ2 dy d|t|) [GeV^-6] of BS09 Eq. 6 for given kinematics (no phase-space checks beyond
  // Q2min); T_pi enters through the elastic cross section.
  double DSigma(int nuPdg, double enu, double Q2, double nu, double absT, double Tpi, int A) const;

  const PionNucleusElastic& Elastic() const { return elastic_; }
  const Options& Opts() const { return opt_; }

 private:
  PionNucleusElastic elastic_;
  Options opt_;
};

}  // namespace nuflow
