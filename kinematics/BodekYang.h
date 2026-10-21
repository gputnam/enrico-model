#pragma once
// Inelastic charged-current structure functions of the nucleon from effective leading-order PDFs:
// the Bodek-Yang model, A. Bodek and U.K. Yang, arXiv:1011.6592 (v2, 2013): GRV98 LO PDFs evaluated
// at the scaling variable xi_w, multiplied by low-Q2 K factors (second-iteration vector parameters of
// Table 2, "type II" axial K factors of Sec. 11.2), with the NMC d/u correction (Sec. 10), the
// H(x) correction of xF3 (Eq. 27), R = sigma_L/sigma_T from R1998 (K. Abe et al., Phys. Lett. B 452,
// 194 (1999)) with the low-Q2 continuation of Sec. 7, and charm production by slow rescaling (Sec. 8).
// See theory/dis.md.

#include <string>

namespace nuflow {

// Dimensionless structure functions F1..F5 of the weak hadronic tensor (Albright-Jarlskog: F4 = 0,
// 2x F5 = F2).
struct StructureFunctions {
  double F1 = 0, F2 = 0, F3 = 0, F4 = 0, F5 = 0;
};

class BodekYang {
 public:
  struct Parameters {
    // xi_w (Eq. 12) [GeV^2]
    double A = 0.621, B = 0.380;
    // vector valence K factors (Eq. 21) [GeV^2]
    double Cv1u = 0.417, Cv2u = 0.264, Cv1d = 0.341, Cv2d = 0.323;
    // vector sea K factors [GeV^2]
    double CseaU = 0.369, CseaD = 0.561, CseaS = 0.561;
    // low-nu factor K_LW = (nu^2 + C)/nu^2 on the vector valence part (W > 1.4 GeV) [GeV^2]
    double ClowNu = 0.218;
    // overall PDF normalisation
    double N = 1.026;
    // type-II axial K factors (Sec. 11.2): (Q2 + P)/(Q2 + C) [GeV^2]
    double PvalAxial = 0.018, CvalAxial = 0.18, PseaAxial = 0.018, CseaAxial = 0.3;
    double MDipoleV2 = 0.71;  // G_D = (1 + Q2/0.71)^-2 in the valence K factor [GeV^2]
    double charmMass = 1.32;  // M_c in xi_w for charm production [GeV]
    bool duCorrection = true;
    bool charm = true;
    bool axialTypeII = true;  // false: type I (axial = vector K factors)
  };

  BodekYang() = default;
  explicit BodekYang(Parameters p) : par_(p) {}

  // nu: true for neutrino (W+ absorbed), false for antineutrino. proton: struck nucleon species.
  // x, Q2 (GeV^2) and nu (energy transfer in the nucleon rest frame, GeV) are those of the hadronic
  // vertex; M is the nucleon mass used in the target-mass terms.
  StructureFunctions Compute(bool neutrino, bool proton, double x, double Q2, double nu, double M) const;

  // xi_w of Eq. 12 for a final-quark mass mf.
  double XiW(double x, double Q2, double M, double mf) const;
  // R = sigma_L / sigma_T used for 2xF1 (R1998 above Q2 = 0.3 GeV^2, continued to 0 at Q2 = 0).
  static double R(double x, double Q2);
  static double R1998(double x, double Q2);
  // H(x) = D_xF3 / D_F2 (Eq. 27).
  static double H(double x) { return 0.914 + 0.296 * x - 0.374 * x * x + 0.165 * x * x * x; }

  const Parameters& Params() const { return par_; }
  std::string Name() const;

 private:
  Parameters par_;
};

}  // namespace nuflow
