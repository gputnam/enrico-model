#pragma once
// GRV98 leading-order parton densities (M. Glueck, E. Reya, A. Vogt, Eur. Phys. J. C 5, 461 (1998),
// hep-ph/9806404), from the authors' interpolation grid (data/pdf/grv98lo.dat). The grid is
// interpolated linearly in (ln x, ln Q2), as in the authors' GRV98PA routine, after dividing out the
// smooth x behaviour of each distribution. Valid for 1e-9 < x < 1 and 0.8 < Q2 < 1e6 GeV^2; Q2 is
// clamped to that range (the Bodek-Yang model freezes the PDFs at Q2 = 0.8 GeV^2).
// See theory/dis.md.

#include <string>
#include <vector>

namespace nuflow {

// x times the distributions in the proton.
struct PartonDensities {
  double uv = 0, dv = 0;      // valence u - ubar, d - dbar
  double ubar = 0, dbar = 0;  // sea antiquarks (= sea quarks of the same flavour)
  double s = 0;               // s = sbar
  double g = 0;               // gluon
};

class GRV98LO {
 public:
  static const GRV98LO& Instance();  // loads data/pdf/grv98lo.dat once
  explicit GRV98LO(const std::string& path);

  PartonDensities Evaluate(double x, double Q2) const;

  static constexpr double kQ2Min = 0.8, kQ2Max = 1e6, kXMin = 1e-9;

 private:
  static constexpr int kNX = 68, kNQ = 27;
  double Interpolate(const std::vector<double>& table, int ix, double dx, int iq, double dq) const;
  std::vector<double> lnx_, lnq2_;
  std::vector<double> uv_, dv_, del_, ud_, s_, g_;
};

}  // namespace nuflow
