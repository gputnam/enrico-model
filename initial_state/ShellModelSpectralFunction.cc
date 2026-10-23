#include "initial_state/ShellModelSpectralFunction.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "core/Constants.h"
#include "core/PDG.h"
#include "fsi/NuclearDensity.h"

namespace nuflow {

using constants::kPi;

namespace {

constexpr double kHbarC = 0.1973269804;  // GeV fm
constexpr int kTable = 20000;

double Gauss(Random& rng) {
  const double u1 = rng.UniformOpen(), u2 = rng.Uniform();
  return std::sqrt(-2 * std::log(u1)) * std::cos(2 * kPi * u2);
}

ThreeVector IsotropicDir(Random& rng) {
  const double c = 2 * rng.Uniform() - 1, s = std::sqrt(std::max(0.0, 1 - c * c));
  const double phi = 2 * kPi * rng.Uniform();
  return {s * std::cos(phi), s * std::sin(phi), c};
}

// |x| for a vector of n independent N(0, sigma^2) components (chi distribution with n dof).
double ChiSample(int n, double sigma, Random& rng) {
  double s = 0;
  for (int i = 0; i < n; ++i) {
    const double g = Gauss(rng);
    s += g * g;
  }
  return sigma * std::sqrt(s);
}

// Ground-state nuclear mass of (Z, A) from the target mass and separation energies [GeV].
// 12C only: residuals after one- and two-nucleon removal (AME2020 mass excesses).
struct TwoNucleonThresholds {
  double Snp = 0.02741;  // 12C -> 10B + n + p
  double Snn = 0.03184;  // 12C -> 10C + n + n
  double Spp = 0.02719;  // 12C -> 10Be + p + p
};

}  // namespace

ShellModelSpectralFunction::ShellModelSpectralFunction() { BuildTables(); }
ShellModelSpectralFunction::ShellModelSpectralFunction(Parameters p) : par_(p) { BuildTables(); }

std::string ShellModelSpectralFunction::Name() const {
  std::ostringstream s;
  s << "SF12C(HO shells, Z_MF=" << par_.mfFraction << ",E1s=" << par_.e1s << ",E1p=" << par_.e1p << ",a=" << par_.widthA
    << ",b=" << par_.widthB << ",Emax=" << par_.eMax << ",kSRC=" << par_.kSRC << "-" << par_.kMax << ",sCM=" << par_.sigmaCM
    << ",np=" << par_.npFraction << (par_.emitPartner ? ",partner" : "") << ")";
  return s.str();
}

void ShellModelSpectralFunction::BuildTables() {
  const double eF = 0.015957;  // 12C proton separation energy: lower edge of the removal-energy spectrum
  const double centroid[2] = {par_.e1s, par_.e1p};
  for (int sh = 0; sh < 2; ++sh) {
    std::vector<double> cdf(kTable + 1, 0.0);
    const double dE = (par_.eMax - eF) / kTable;
    for (int i = 0; i < kTable; ++i) {
      const double E = eF + (i + 0.5) * dE;
      const double x = E - eF;
      const double G = par_.widthA * x * x / (par_.widthB + x * x);
      const double d = E - centroid[sh];
      // Lorentzian with energy-dependent width; a minimum width of 10 keV keeps the 1p peak resolvable.
      const double Ge = std::max(G, 1e-5);
      cdf[i + 1] = cdf[i] + Ge / (2 * kPi) / (d * d + 0.25 * Ge * Ge) * dE;
    }
    for (auto& c : cdf) c /= cdf.back();
    // Inverse CDF on a uniform grid in u.
    eGrid_[sh].assign(kTable + 1, 0.0);
    int j = 0;
    for (int i = 0; i <= kTable; ++i) {
      const double u = static_cast<double>(i) / kTable;
      while (j < kTable && cdf[j + 1] < u) ++j;
      const double c0 = cdf[j], c1 = cdf[std::min(j + 1, kTable)];
      const double f = c1 > c0 ? (u - c0) / (c1 - c0) : 0.0;
      eGrid_[sh][i] = eF + (j + std::clamp(f, 0.0, 1.0)) * dE;
    }
  }
}

double ShellModelSpectralFunction::SampleRemovalEnergy(int shell, double u) const {
  const double x = u * kTable;
  const int i = std::min(static_cast<int>(x), kTable - 1);
  const double f = x - i;
  return eGrid_[shell][i] * (1 - f) + eGrid_[shell][i + 1] * f;
}

void ShellModelSpectralFunction::Sample(const Nucleus& nucleus, const HitNucleonRequest& request, InitialState& out,
                                        Random& rng) const {
  if (nucleus.Z != 6 || nucleus.A != 12)
    throw std::runtime_error("ShellModelSpectralFunction: parameters are only available for 12C");
  if (request.nNucleons != 1) throw std::logic_error("ShellModelSpectralFunction: single-nucleon channels only");
  const NuclearDensity& density = NuclearDensity::For(nucleus);
  const double b = density.B();            // HO length [fm]
  const double sigK = kHbarC / b / std::sqrt(2.0);  // momentum-space Gaussian sigma per component [GeV]
  const double sigR = b / std::sqrt(2.0);           // coordinate-space sigma per component [fm]
  out.nHit = 1;
  auto& h = out.hit[0];
  h.pdg = request.species[0];
  const double mHit = pdg::NucleonMass(h.pdg);
  const double dS = h.pdg == pdg::kNeutron ? nucleus.separationNeutron - nucleus.separationProton : 0.0;
  out.hasPosition = true;
  out.hasSpectator = false;

  if (rng.Uniform() < par_.mfFraction) {
    // Mean field: 1s (2 per species) or 1p (Z - 2 per species).
    const double pS = 2.0 / nucleus.Z;
    const int shell = rng.Uniform() < pS ? 0 : 1;
    const int dof = shell == 0 ? 3 : 5;  // |phi_1s|^2 ~ e^{-y^2}, |phi_1p|^2 ~ y^2 e^{-y^2}
    h.p = IsotropicDir(rng) * ChiSample(dof, sigK, rng);
    out.position = IsotropicDir(rng) * ChiSample(dof, sigR, rng);
    const double ER = SampleRemovalEnergy(shell, rng.Uniform()) + dS;
    const double mRes = nucleus.mass - mHit + ER;  // excited A-1 residual
    h.energy = nucleus.mass - std::sqrt(mRes * mRes + h.p.Mag2());
    h.removalEnergy = ER;
    FillRemnant(nucleus, out);
    return;
  }

  // Short-range-correlated pair.
  out.position = density.SamplePosition(rng);
  const double u = rng.Uniform();
  const double k = 1.0 / (1.0 / par_.kSRC - u * (1.0 / par_.kSRC - 1.0 / par_.kMax));  // k^2 n(k) ~ k^-2
  h.p = IsotropicDir(rng) * k;
  const bool other = rng.Uniform() < par_.npFraction;
  const int partner = other ? (h.pdg == pdg::kProton ? pdg::kNeutron : pdg::kProton) : h.pdg;
  const ThreeVector Kcm{par_.sigmaCM * Gauss(rng), par_.sigmaCM * Gauss(rng), par_.sigmaCM * Gauss(rng)};
  const ThreeVector p2 = Kcm - h.p;
  const TwoNucleonThresholds thr;
  const double S2 = other ? thr.Snp : (h.pdg == pdg::kNeutron ? thr.Snn : thr.Spp);
  const double m2 = pdg::NucleonMass(partner);
  const double mRes2 = nucleus.mass - mHit - m2 + S2;  // A-2 residual in its ground state
  const double e2 = std::sqrt(m2 * m2 + p2.Mag2());
  const double eRes = std::sqrt(mRes2 * mRes2 + Kcm.Mag2());
  h.energy = nucleus.mass - e2 - eRes;
  h.removalEnergy = mHit - h.energy;  // includes the partner's kinetic energy and the pair recoil
  if (par_.emitPartner) {
    out.hasSpectator = true;
    out.spectator.pdg = partner;
    out.spectator.p = p2;
    out.spectator.energy = e2;
    out.spectator.removalEnergy = 0;
  } else {
    // The partner stays bound: equivalent removal from an excited A-1 residual with the same E_R.
    const double mRes = nucleus.mass - mHit + h.removalEnergy;
    h.energy = nucleus.mass - std::sqrt(mRes * mRes + h.p.Mag2());
  }
  FillRemnant(nucleus, out);
}

bool ShellModelSpectralFunction::PauliBlocked(const Nucleus&, int, const ThreeVector&) const {
  throw std::logic_error("ShellModelSpectralFunction: Pauli blocking needs the interaction point");
}

bool ShellModelSpectralFunction::PauliBlocked(const Nucleus& nucleus, int pdgCode, const ThreeVector& p,
                                              const InitialState& initial) const {
  const double kF = NuclearDensity::For(nucleus).KFermi(pdgCode, initial.position.Mag());
  return p.Mag2() < kF * kF;
}

MediumPoint ShellModelSpectralFunction::Medium(const Nucleus& nucleus, const InitialState& initial) const {
  const NuclearDensity& density = NuclearDensity::For(nucleus);
  MediumPoint m;
  m.rho = density.Rho(initial.position.Mag());
  m.rhoCentral = density.RhoCentral();
  return m;
}

}  // namespace nuflow
