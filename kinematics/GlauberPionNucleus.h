#pragma once
// Glauber (eikonal multiple-scattering) pion-nucleus elastic amplitude, used by the coherent pion channel
// above T_pi ~ 0.9 GeV (theory/coherent_pion.md, section 1.1).
//
//   F(q)   = i k Int_0^inf b db J0(q b) Gamma(b) * exp(q^2 b_HO^2 / 4A)            (c.m. correction, A <= 16)
//   Gamma(b) = 1 - [1 - sigma_piN (1 - i rho) / 2 * T_beta(b)]^A                 (Glauber 1959; Czyz-Maximon 1969)
//   T_beta(b) = (1/2pi) Int_0^inf q dq J0(q b) S(q) exp(-beta q^2 / 2)            (thickness folded with the piN range)
//   dsigma/dt = pi |Int b db J0(q b) Gamma(b)|^2,  q = sqrt|t|
//
// Inputs, none fitted to neutrino data:
//  - S(q): point-nucleon single-particle form factor of the ground state. 12C and 16O: harmonic-oscillator
//    charge densities from elastic electron scattering (H. de Vries, C.W. de Jager, C. de Vries, At. Data
//    Nucl. Data Tables 36, 495 (1987)), the nucleon charge radii (r_p^2 + r_n^2) unfolded, and the HO
//    fixed-centre (c.m.) factor. Heavier nuclei: two-parameter Fermi (de Vries 1987), used as point density.
//  - sigma_piN(p): isoscalar (pi+ p + pi- p)/2 total cross section from the PDG data compilation (RPP 2022
//    files), including the N* region, Fermi-averaged in s with the 12C Fermi momentum.
//  - rho(p) = Re f / Im f (forward, isoscalar): derivative-analyticity of the PDG Regge + ln^2 s fit.
//  - beta(p): piN slope from the measured elastic cross section, sigma_el = sigma_tot^2 (1 + rho^2)/(16 pi beta).

#include <array>
#include <atomic>
#include <complex>
#include <memory>
#include <mutex>
#include <vector>

namespace nuflow {

class GlauberPionNucleus {
 public:
  struct PiN {
    double sigmaTot = 0;  // isoscalar total [mb]
    double sigmaEl = 0;   // isoscalar elastic [mb]
    double rho = 0;       // Re/Im of the forward amplitude
    double beta = 0;      // slope of the amplitude, f ~ exp(-beta q^2 / 2) [GeV^-2]
  };
  // Isoscalar pi N inputs at pion lab momentum p [GeV].
  static PiN PiNAmplitude(double p);
  // rho from the PDG Regge fit (derivative analyticity), exposed for validation.
  static double ReggeRho(double p);

  struct Observables {
    double sigmaTot = 0;   // pi A total [mb]
    double sigmaEl = 0;    // coherent elastic [mb]
    double sigmaInel = 0;  // sigmaTot - sigmaEl (includes quasi-elastic) [mb]
    double sigmaProd = 0;  // particle production ("absorption"): Int d2b {1 - [1 - sigma_in T]^A} [mb]
  };

  GlauberPionNucleus();
  // dsigma/dt [mb/GeV^2] at pion kinetic energy T [GeV] and |t| [GeV^2] (interpolated on a precomputed
  // grid; built on first use for each A; thread safe).
  double DSigmaDt(double T, double absT, int A) const;
  // Direct (non-tabulated) evaluation, for validation.
  Observables Compute(double T, int A) const;
  double DSigmaDtDirect(double T, double absT, int A) const;

  static constexpr double kTMin = 0.6, kTMax = 250.0;  // tabulated T range [GeV]

 private:
  struct Density;
  struct Table {
    std::vector<double> lnT;                      // grid in ln T
    std::vector<std::vector<std::complex<double>>> amp;  // [iT][iq] amplitude Int b J0 Gamma db [fm^2]
  };
  const Table& GetTable(int A) const;
  std::vector<double> Thickness(double T, int A, PiN& pin) const;  // T_beta(b) [fm^-2], normalised to 1
  std::vector<std::complex<double>> Profile(double T, int A, PiN* used = nullptr) const;  // Gamma(b)

  // b grid [fm] and q grid of the amplitude table (q = sqrt|t|) [GeV]
  static constexpr int kNb = 800;
  static constexpr double kBMax = 16.0;
  static constexpr int kNq = 321;
  static constexpr double kQMax = 1.0;  // |t| up to 1 GeV^2
  static constexpr int kNT = 200;

  struct Cache;
  std::shared_ptr<Cache> cache_;
};

}  // namespace nuflow
