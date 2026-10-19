#pragma once
// Nucleon form factors for the charged-current (isovector) weak current.
//
// Convention: F_A(0) = +g_A > 0 with the hadronic current
//   J^mu = ubar(p') [ gamma^mu F1 + i sigma^{mu nu} q_nu F2 / 2M + gamma^mu gamma5 F_A + q^mu gamma5 F_P / M ] u(p),
// the neutrino/antineutrino difference entering through the lepton tensor (see theory/quasielastic.md).
// F2 here already includes the anomalous moment (F2(0) = mu_p - mu_n - 1 = 3.706).

#include <array>
#include <memory>
#include <string>

#include <yaml-cpp/yaml.h>

namespace nuflow {

struct SachsFormFactors {
  double GEp = 0, GMp = 0, GEn = 0, GMn = 0;
};

struct WeakFormFactors {
  double F1 = 0, F2 = 0, FA = 0, FP = 0;
};

// Electromagnetic Sachs form factors of the proton and neutron.
class VectorFormFactorModel {
 public:
  virtual ~VectorFormFactorModel() = default;
  virtual std::string Name() const = 0;
  virtual SachsFormFactors Sachs(double Q2) const = 0;
};

// Dipole G_Ep, G_Mp, G_Mn with Galster G_En.
class DipoleGalsterFormFactors : public VectorFormFactorModel {
 public:
  explicit DipoleGalsterFormFactors(double MV2 = 0.71) : MV2_(MV2) {}
  std::string Name() const override { return "dipole+Galster"; }
  SachsFormFactors Sachs(double Q2) const override;

 private:
  double MV2_;
};

// Bradford, Bodek, Budd & Arrington, Nucl. Phys. B Proc. Suppl. 159 (2006) 127 [hep-ex/0602017].
class BBBA05FormFactors : public VectorFormFactorModel {
 public:
  std::string Name() const override { return "BBBA05"; }
  SachsFormFactors Sachs(double Q2) const override;
};

class AxialFormFactorModel {
 public:
  virtual ~AxialFormFactorModel() = default;
  virtual std::string Name() const = 0;
  virtual double FA(double Q2) const = 0;
};

class DipoleAxialFormFactor : public AxialFormFactorModel {
 public:
  DipoleAxialFormFactor(double gA, double MA) : gA_(gA), MA2_(MA * MA) {}
  std::string Name() const override { return "dipole"; }
  double FA(double Q2) const override {
    const double d = 1.0 + Q2 / MA2_;
    return gA_ / (d * d);
  }

 private:
  double gA_, MA2_;
};

// Model-independent z expansion: F_A = sum_k a_k z^k with
//   z = (sqrt(tcut + Q2) - sqrt(tcut - t0)) / (sqrt(tcut + Q2) + sqrt(tcut - t0)).
// The free parameters a_1..a_4 are given; a_0 fixes F_A(0) = g_A and a_5..a_8 are fixed by the
// four sum rules sum_k k^n a_k = 0 (n = 0..3) that impose F_A ~ 1/Q^4 at large Q2
// (Meyer, Betancourt, Gran & Hill, Phys. Rev. D 93, 113015 (2016) [arXiv:1603.03048]).
class ZExpansionAxialFormFactor : public AxialFormFactorModel {
 public:
  ZExpansionAxialFormFactor(double gA, std::array<double, 4> a1to4, double tcut, double t0);
  std::string Name() const override { return "z-expansion"; }
  double FA(double Q2) const override;
  const std::array<double, 9>& Coefficients() const { return a_; }

 private:
  double Z(double Q2) const;
  std::array<double, 9> a_{};
  double tcut_, t0_;
};

// Complete set of CC form factors: CVC isovector vector part + axial + PCAC pseudoscalar.
class NucleonFormFactors {
 public:
  NucleonFormFactors(std::shared_ptr<const VectorFormFactorModel> vec, std::shared_ptr<const AxialFormFactorModel> ax)
      : vector_(std::move(vec)), axial_(std::move(ax)) {}
  WeakFormFactors Evaluate(double Q2) const;
  std::string Name() const { return vector_->Name() + "/" + axial_->Name(); }

  static NucleonFormFactors FromConfig(const YAML::Node& node);

 private:
  std::shared_ptr<const VectorFormFactorModel> vector_;
  std::shared_ptr<const AxialFormFactorModel> axial_;
};

}  // namespace nuflow
