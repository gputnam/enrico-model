#include "fsi/IntranuclearCascade.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <memory>

#include "core/Constants.h"
#include "core/PDG.h"
#include "fsi/HadronNucleonXSec.h"

namespace nuflow {

using constants::kPi;

namespace {

constexpr double kHbarC = 0.1973269804;    // GeV fm
constexpr double kGeV2mb = 0.389379;       // 1 GeV^-2 = 0.389379 mb
constexpr double kMbToFm2 = 0.1;

bool IsNucleon(int c) { return c == pdg::kProton || c == pdg::kNeutron; }
bool IsPion(int c) { return c == pdg::kPiPlus || c == pdg::kPiMinus || c == pdg::kPiZero; }
int Charge(int c) {
  switch (c) {
    case pdg::kProton: case pdg::kPiPlus: return 1;
    case pdg::kPiMinus: return -1;
    default: return 0;
  }
}
double Mass(int c) { return IsNucleon(c) ? pdg::NucleonMass(c) : pdg::PionMass(c); }
int PionWithCharge(int q) { return q > 0 ? pdg::kPiPlus : (q < 0 ? pdg::kPiMinus : pdg::kPiZero); }
int NucleonWithCharge(int q) { return q > 0 ? pdg::kProton : pdg::kNeutron; }

ThreeVector IsotropicDirection(Random& rng) {
  const double c = 2 * rng.Uniform() - 1, s = std::sqrt(std::max(0.0, 1 - c * c));
  const double phi = 2 * kPi * rng.Uniform();
  return {s * std::cos(phi), s * std::sin(phi), c};
}

// Unit vector at polar angle acos(cosT) and azimuth phi around `axis`.
ThreeVector AroundAxis(const ThreeVector& axis, double cosT, double phi) {
  const ThreeVector z = axis.Unit();
  const ThreeVector a = std::abs(z.x) < 0.9 ? ThreeVector{1, 0, 0} : ThreeVector{0, 1, 0};
  const ThreeVector x = z.Cross(a).Unit();
  const ThreeVector y = z.Cross(x);
  const double sinT = std::sqrt(std::max(0.0, 1 - cosT * cosT));
  return z * cosT + x * (sinT * std::cos(phi)) + y * (sinT * std::sin(phi));
}

// cos(theta) with density exp(a cos(theta)) on [-1, 1].
double SampleExpCos(double a, Random& rng) {
  const double u = rng.Uniform();
  if (a < 1e-6) return 2 * u - 1;
  return std::clamp(1 + std::log(u + (1 - u) * std::exp(-2 * a)) / a, -1.0, 1.0);
}

FourVector OnShell(double m, const ThreeVector& p) { return {std::sqrt(p.Mag2() + m * m), p}; }

// Fermi-sea nucleon: |k| uniform in the sphere of radius kF.
FourVector FermiSeaNucleon(int code, double kF, Random& rng) {
  const double k = kF * std::cbrt(rng.UniformOpen());
  return OnShell(pdg::NucleonMass(code), IsotropicDirection(rng) * k);
}

// Two-body final state of total momentum P: masses m3, m4, particle 3 along CM direction n.
void TwoBody(const FourVector& P, double m3, double m4, const ThreeVector& nCM, FourVector& p3, FourVector& p4) {
  const double W = P.M();
  const double q = hadxs::CMMomentum(W, m3, m4);
  const ThreeVector b = P.BoostVector();
  p3 = OnShell(m3, nCM * q).Boosted(b);
  p4 = OnShell(m4, nCM * (-q)).Boosted(b);
}

// Moller velocity of the pair divided by the projectile velocity (rate per unit path length).
double FluxFactor(const FourVector& p1, const FourVector& p2) {
  const double m1 = p1.M(), m2 = p2.M();
  const double d = p1.Dot(p2);
  const double vMol = std::sqrt(std::max(0.0, d * d - m1 * m1 * m2 * m2)) / (p1.e * p2.e);
  const double v1 = p1.P() / p1.e;
  return v1 > 1e-9 ? vMol / v1 : 0.0;
}

// Momentum of particle 1 in the rest frame of particle 2.
double LabMomentum(const FourVector& p1, const FourVector& p2) {
  const double m1 = p1.M(), m2 = p2.M();
  const double e = p1.Dot(p2) / m2;
  return std::sqrt(std::max(0.0, e * e - m1 * m1));
}

// Fraction of isotropic (Delta rest frame) decays Delta -> N pi whose nucleon ends up above the Fermi
// momentum kF in the frame where the Delta has four-momentum P: the Pauli-blocked width is
// Gamma~ = Gamma * PauliFraction (Oset & Salcedo use the same Gamma~ in the Delta propagator).
double PauliFraction(const FourVector& P, double mN, double mPi, double kF) {
  const double W = P.M();
  const double q = hadxs::CMMomentum(W, mN, mPi);
  if (q <= 0) return 0.0;
  const double eN = std::sqrt(q * q + mN * mN);
  const double gamma = P.e / W, gb = P.P() / W;  // gamma*beta
  // |p_lab|^2 = (gamma*beta)^2 q^2 c^2 + 2 gamma^2 beta eN q c + (gamma*beta)^2 eN^2 + q^2, increasing in c
  auto p2 = [&](double c) { return gb * gb * q * q * c * c + 2 * gamma * gb * eN * q * c + gb * gb * eN * eN + q * q; };
  const double k2 = kF * kF;
  if (p2(-1) >= k2) return 1.0;
  if (p2(1) <= k2) return 0.0;
  double lo = -1, hi = 1;
  for (int i = 0; i < 40; ++i) {
    const double mid = 0.5 * (lo + hi);
    (p2(mid) < k2 ? lo : hi) = mid;
  }
  return 0.5 * (1 - 0.5 * (lo + hi));
}

}  // namespace

const NuclearDensity& IntranuclearCascade::DensityFor(const Nucleus& nucleus) {
  thread_local std::map<int, std::unique_ptr<NuclearDensity>> cache;
  auto it = cache.find(nucleus.pdg);
  if (it == cache.end()) it = cache.emplace(nucleus.pdg, std::make_unique<NuclearDensity>(nucleus)).first;
  return *it->second;
}

double IntranuclearCascade::Potential(const Nucleus& nucleus, const NuclearDensity& density, int code,
                                      double r) const {
  const double m = pdg::NucleonMass(code);
  const double kF = density.KFermi(code, r);
  double S = code == pdg::kProton ? nucleus.separationProton : nucleus.separationNeutron;
  if (S <= 0) S = 0.008;
  return std::sqrt(kF * kF + m * m) - m + S;
}

void IntranuclearCascade::AddNucleon(const Nucleus& nucleus, const NuclearDensity& density, int code,
                                     const ThreeVector& x, const FourVector& p, int ancestor, int generation,
                                     std::vector<Hadron>& active, Remnant& rem, Stats* stats) const {
  const double T = p.e - pdg::NucleonMass(code);
  const double estar = T - Potential(nucleus, density, code, x.Mag());
  if (estar <= 0) {  // bound: joins the remnant
    if (code == pdg::kProton) ++rem.z;
    else ++rem.n;
    if (stats) ++stats->captured;
    return;
  }
  Hadron h;
  h.pdg = code;
  h.x = x;
  h.dir = p.Vect().Unit();
  h.estar = estar;
  h.ancestor = ancestor;
  h.generation = generation;
  active.push_back(h);
}

void IntranuclearCascade::Propagate(const Nucleus& nucleus, Event& event, Random& rng) const {
  auto& parts = event.particles;
  const size_t n = parts.size();
  if (nucleus.IsFreeNucleon()) {  // hydrogen: nothing to rescatter on
    for (size_t i = 0; i < n; ++i) {
      const Particle part = parts[i];
      if (part.status == kStatusPreFSI) event.Add(part.pdg, kStatusFinal, part.p, static_cast<int>(i));
    }
    return;
  }
  const NuclearDensity& density = DensityFor(nucleus);
  // The primary vertex: from the stage-2 nuclear model if it has a radial dependence (local Fermi gas),
  // otherwise sampled from the density.
  const ThreeVector vertex = event.initial.hasPosition ? event.initial.position : density.SamplePosition(rng);
  const bool formation =
      opt_.formationZone && std::find(opt_.formationProcesses.begin(), opt_.formationProcesses.end(),
                                      event.channel) != opt_.formationProcesses.end();

  int zRem = nucleus.Z, nRem = nucleus.N();
  std::vector<Hadron> active;
  std::vector<int> passThrough;
  for (size_t i = 0; i < n; ++i) {
    Particle& part = parts[i];
    if (part.status == kStatusRemnant) {
      zRem = pdg::IonZ(part.pdg);
      nRem = pdg::IonA(part.pdg) - zRem;
    }
    if (part.status != kStatusPreFSI) continue;
    part.x = vertex;
    if (!IsNucleon(part.pdg) && !IsPion(part.pdg)) {
      passThrough.push_back(static_cast<int>(i));
      continue;
    }
    Hadron h;
    h.pdg = part.pdg;
    h.x = vertex;
    h.dir = part.p.Vect().Unit();
    h.estar = std::max(0.0, part.p.e - Mass(part.pdg));
    h.ancestor = static_cast<int>(i);
    if (formation) h.formation = part.p.P() / opt_.formationMu2 * kHbarC;
    active.push_back(h);
  }
  std::vector<Hadron> escaped;
  Cascade(nucleus, density, std::move(active), escaped, zRem, nRem, rng);

  // Final state: escaping hadrons, then the remnant with the balancing four-momentum.
  FourVector in, out;
  for (size_t i = 0; i < n; ++i) {
    if (parts[i].status == kStatusBeam || parts[i].status == kStatusTarget) in = in + parts[i].p;
    if (parts[i].status == kStatusFinal) out = out + parts[i].p;
  }
  for (int i : passThrough) {
    const Particle part = parts[i];
    const int k = event.Add(part.pdg, kStatusFinal, part.p, i);
    parts[k].x = vertex;
    out = out + part.p;
  }
  for (const auto& h : escaped) {
    const double m = Mass(h.pdg);
    const double p = std::sqrt(h.estar * (h.estar + 2 * m));
    const FourVector p4(m + h.estar, h.dir * p);
    const int k = event.Add(h.pdg, kStatusFinal, p4, h.ancestor);
    parts[k].x = h.x;
    out = out + p4;
  }
  const int aRem = zRem + nRem;
  if (aRem > 0) event.Add(pdg::IonCode(zRem, aRem), kStatusRemnantPostFSI, in - out, -1);
}

void IntranuclearCascade::Cascade(const Nucleus& nucleus, const NuclearDensity& density, std::vector<Hadron> active,
                                  std::vector<Hadron>& escaped, int& zRem, int& nRem, Random& rng,
                                  Stats* stats) const {
  Remnant rem{zRem, nRem, nucleus.Z, nucleus.N()};
  int guard = 0;
  while (!active.empty()) {
    Hadron h = active.back();
    active.pop_back();
    if (++guard > 10000) {  // safety net, never reached in practice
      escaped.push_back(h);
      continue;
    }
    Transport(nucleus, density, h, active, escaped, rem, rng, stats);
  }
  zRem = rem.z;
  nRem = rem.n;
}

void IntranuclearCascade::Transport(const Nucleus& nucleus, const NuclearDensity& density, Hadron h,
                                    std::vector<Hadron>& active, std::vector<Hadron>& escaped, Remnant& rem,
                                    Random& rng, Stats* stats) const {
  const double R = density.RMax();
  const bool nucleon = IsNucleon(h.pdg);
  for (int istep = 0; istep < 100000; ++istep) {
    const double r = h.x.Mag();
    if (r >= R && h.x.Dot(h.dir) >= 0) {
      escaped.push_back(h);
      return;
    }
    double ds = opt_.step;
    if (h.formation > 0) {
      ds = std::min(ds, h.formation);
      h.formation -= ds;
      h.x = h.x + h.dir * ds;
      continue;
    }
    bool consumed = false;
    const bool interacted = nucleon ? NucleonInteraction(nucleus, density, h, r, active, rem, rng, stats, consumed)
                                    : PionInteraction(nucleus, density, h, r, active, rem, rng, stats, consumed);
    if (interacted && consumed) return;  // products (if any) are in `active`
    h.x = h.x + h.dir * ds;
  }
  escaped.push_back(h);
}

// Nucleon: decides whether a collision happens in the next step and performs it.
bool IntranuclearCascade::NucleonInteraction(const Nucleus& nucleus, const NuclearDensity& density, Hadron& h,
                                             double r, std::vector<Hadron>& active, Remnant& rem, Random& rng,
                                             Stats* stats, bool& consumed) const {
  consumed = false;
  if (!opt_.nucleonRescattering) return false;
  const double m1 = pdg::NucleonMass(h.pdg);
  const double tIn = h.estar + Potential(nucleus, density, h.pdg, r);
  const FourVector p1(m1 + tIn, h.dir * std::sqrt(tIn * (tIn + 2 * m1)));
  FourVector p2[2];
  double rate[2] = {0, 0}, sigTot[2] = {0, 0};
  const int species[2] = {pdg::kProton, pdg::kNeutron};
  // A knock-out must leave at least one nucleon in the remnant (it carries the binding energy).
  if (rem.z + rem.n < 2) return false;
  for (int t = 0; t < 2; ++t) {
    const double frac = rem.Frac(t == 0);
    if (frac <= 0) continue;
    p2[t] = FermiSeaNucleon(species[t], density.KFermi(species[t], r), rng);
    const bool like = species[t] == h.pdg;
    sigTot[t] = hadxs::NNTotal(like, LabMomentum(p1, p2[t]));
    rate[t] = density.RhoSpecies(species[t], r) * frac * sigTot[t] * kMbToFm2 * FluxFactor(p1, p2[t]);
  }
  const double lam = rate[0] + rate[1];
  if (lam <= 0 || rng.Uniform() >= 1 - std::exp(-lam * opt_.step)) return false;
  const int t = rng.Uniform() * lam < rate[0] ? 0 : 1;
  const int code2 = species[t];
  const bool like = code2 == h.pdg;
  const FourVector P = p1 + p2[t];
  const double W = P.M();
  const double plab = LabMomentum(p1, p2[t]);
  const double sigEl = hadxs::NNElastic(like, plab);
  const ThreeVector b = P.BoostVector();
  const ThreeVector axis = p1.Boosted(-b).Vect();
  const double mPi = constants::kMassPionCharged;
  const bool inelastic = opt_.nnPionProduction && rng.Uniform() * sigTot[t] >= sigEl &&
                         W > 2 * constants::kMassNeutron + mPi + 0.002;
  if (stats) ++stats->nucleonCollisions;

  if (!inelastic) {
    const double q = hadxs::CMMomentum(W, m1, pdg::NucleonMass(code2));
    double c = SampleExpCos(2 * hadxs::NNSlope(like, plab) * q * q, rng);
    if (rng.Uniform() < 0.5) c = -c;  // forward/backward symmetric (identical particles; np charge exchange)
    const ThreeVector n = AroundAxis(axis, c, 2 * kPi * rng.Uniform());
    FourVector p3, p4;
    TwoBody(P, m1, pdg::NucleonMass(code2), n, p3, p4);
    if (p3.P() < density.KFermi(h.pdg, r) || p4.P() < density.KFermi(code2, r)) {
      if (stats) ++stats->nucleonPauliBlocked;
      return false;
    }
    if (code2 == pdg::kProton) --rem.z;
    else --rem.n;
    AddNucleon(nucleus, density, h.pdg, h.x, p3, h.ancestor, h.generation + 1, active, rem, stats);
    AddNucleon(nucleus, density, code2, h.x, p4, h.ancestor, h.generation + 1, active, rem, stats);
    consumed = true;
    return true;
  }

  // N N -> N Delta (isospin Clebsch-Gordan), Delta -> N pi.
  const int qTot = Charge(h.pdg) + Charge(code2);
  int qN, qD;
  if (qTot == 2) {        // pp: n Delta++ (3/4), p Delta+ (1/4)
    if (rng.Uniform() < 0.75) qN = 0, qD = 2;
    else qN = 1, qD = 1;
  } else if (qTot == 0) { // nn: p Delta- (3/4), n Delta0 (1/4)
    if (rng.Uniform() < 0.75) qN = 1, qD = -1;
    else qN = 0, qD = 0;
  } else {                // np: p Delta0 (1/2), n Delta+ (1/2)
    if (rng.Uniform() < 0.5) qN = 1, qD = 0;
    else qN = 0, qD = 1;
  }
  int qNd, qPi;  // Delta decay charges
  switch (qD) {
    case 2: qNd = 1, qPi = 1; break;
    case 1: if (rng.Uniform() < 2.0 / 3.0) qNd = 1, qPi = 0; else qNd = 0, qPi = 1; break;
    case 0: if (rng.Uniform() < 2.0 / 3.0) qNd = 0, qPi = 0; else qNd = 1, qPi = -1; break;
    default: qNd = 0, qPi = -1; break;
  }
  const int cN = NucleonWithCharge(qN), cNd = NucleonWithCharge(qNd), cPi = PionWithCharge(qPi);
  const double mN = pdg::NucleonMass(cN), mNd = pdg::NucleonMass(cNd), mP = pdg::PionMass(cPi);
  const double mMin = mNd + mP + 1e-4, mMax = W - mN - 1e-4;
  if (mMax <= mMin) return false;
  // Delta mass: relativistic Breit-Wigner (energy-dependent width) x two-body phase space.
  auto density2 = [&](double m) {
    const double G = hadxs::DeltaWidth(m), MD = hadxs::kMassDelta;
    return m * m * G / ((m * m - MD * MD) * (m * m - MD * MD) + MD * MD * G * G) * hadxs::CMMomentum(W, mN, m);
  };
  double fmax = 0;
  for (int i = 0; i <= 40; ++i) fmax = std::max(fmax, density2(mMin + (mMax - mMin) * i / 40.0));
  fmax *= 1.3;
  double mD = 0;
  for (int it = 0; it < 1000; ++it) {
    const double m = rng.Uniform(mMin, mMax);
    if (rng.Uniform() * fmax < density2(m)) {
      mD = m;
      break;
    }
  }
  if (mD <= 0) return false;
  const double qi = hadxs::CMMomentum(W, m1, pdg::NucleonMass(code2)), qf = hadxs::CMMomentum(W, mN, mD);
  double c = SampleExpCos(2 * hadxs::NNSlope(like, plab) * qi * qf, rng);
  if (rng.Uniform() < 0.5) c = -c;
  FourVector pN, pD;
  TwoBody(P, mN, mD, AroundAxis(axis, c, 2 * kPi * rng.Uniform()), pN, pD);
  FourVector pNd, pPi;
  TwoBody(pD, mNd, mP, IsotropicDirection(rng), pNd, pPi);
  if (pN.P() < density.KFermi(cN, r) || pNd.P() < density.KFermi(cNd, r)) {
    if (stats) ++stats->nucleonPauliBlocked;
    return false;
  }
  if (code2 == pdg::kProton) --rem.z;
  else --rem.n;
  if (stats) ++stats->nnInelastic;
  AddNucleon(nucleus, density, cN, h.x, pN, h.ancestor, h.generation + 1, active, rem, stats);
  AddNucleon(nucleus, density, cNd, h.x, pNd, h.ancestor, h.generation + 1, active, rem, stats);
  Hadron pi;
  pi.pdg = cPi;
  pi.x = h.x;
  pi.dir = pPi.Vect().Unit();
  pi.estar = pPi.e - mP;
  pi.ancestor = h.ancestor;
  pi.generation = h.generation + 1;
  active.push_back(pi);
  consumed = true;
  return true;
}

// Pion: Delta-hole model with the Oset-Salcedo in-medium Delta width.
bool IntranuclearCascade::PionInteraction(const Nucleus& nucleus, const NuclearDensity& density, Hadron& h, double r,
                                          std::vector<Hadron>& active, Remnant& rem, Random& rng, Stats* stats,
                                          bool& consumed) const {
  consumed = false;
  const double mPi = pdg::PionMass(h.pdg);
  const double tPi = h.estar;
  if (tPi <= 0) return false;
  const FourVector pPi(mPi + tPi, h.dir * std::sqrt(tPi * (tPi + 2 * mPi)));
  const int qPi = Charge(h.pdg);
  const double rho = density.Rho(r);
  const hadxs::OsetSelfEnergy se = hadxs::DeltaSelfEnergy(tPi, rho / hadxs::kRho0);
  const int species[2] = {pdg::kProton, pdg::kNeutron};
  FourVector p2[2];
  // channel rates per species: 0 QE (Delta decay + C_Q), 1 non-resonant, 2 abs2, 3 abs3
  double rate[2][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}};
  double lam = 0;
  for (int t = 0; t < 2; ++t) {
    const double frac = rem.Frac(t == 0);
    if (frac <= 0) continue;
    p2[t] = FermiSeaNucleon(species[t], density.KFermi(species[t], r), rng);
    const FourVector P = pPi + p2[t];
    const double W = P.M();
    const double q = hadxs::CMMomentum(W, pdg::NucleonMass(species[t]), mPi);
    if (q <= 0) continue;
    const double wD = hadxs::DeltaRegionWeight(W);
    double sig[4] = {0, 0, 0, 0};
    if (wD > 0) {
      const double G = hadxs::DeltaWidth(W);
      // Pauli-blocked decay width (decay nucleon above the local Fermi sea)
      const double Gt = G * PauliFraction(P, pdg::NucleonMass(species[t]), mPi, density.KFermi(species[t], r));
      const double cA2 = opt_.pionAbsorption ? se.cA2 : 0, cA3 = opt_.pionAbsorption ? se.cA3 : 0;
      const double Gtot = Gt + 2 * (se.cQ + se.cA2 + se.cA3);
      const double D = (W - hadxs::kMassDelta) * (W - hadxs::kMassDelta) + 0.25 * Gtot * Gtot;
      // formation with the free coupling (Gamma/2), decay into channel c with Gamma_c/2
      const double pref = wD * hadxs::DeltaIsospinFactor(qPi, Charge(species[t])) * 8 * kPi / (q * q) * kGeV2mb *
                          (0.5 * G) / D;
      sig[0] = opt_.pionQuasiElastic ? pref * (0.5 * Gt + se.cQ) : 0;
      sig[2] = pref * cA2;
      sig[3] = pref * cA3;
    }
    sig[1] = opt_.pionQuasiElastic ? (1 - wD) * hadxs::kPiNHighEnergy : 0;
    const double f = density.RhoSpecies(species[t], r) * frac * kMbToFm2 * FluxFactor(pPi, p2[t]);
    for (int c = 0; c < 4; ++c) {
      rate[t][c] = f * sig[c];
      lam += rate[t][c];
    }
  }
  if (lam <= 0 || rng.Uniform() >= 1 - std::exp(-lam * opt_.step)) return false;
  double u = rng.Uniform() * lam;
  int t = 0, ch = 0;
  for (t = 0; t < 2; ++t) {
    for (ch = 0; ch < 4; ++ch) {
      if (u < rate[t][ch]) break;
      u -= rate[t][ch];
    }
    if (ch < 4) break;
  }
  if (t >= 2) t = 1, ch = 3;
  const int c1 = species[t];
  const int q1 = Charge(c1);
  const FourVector P = pPi + p2[t];
  const double W = P.M();

  if (rem.z + rem.n < (ch <= 1 ? 2 : ch) + 1) return false;  // keep at least one nucleon in the remnant
  if (ch <= 1) {  // quasi-elastic: Delta decay (or non-resonant elastic-like scattering)
    int qPiOut = qPi, qNOut = q1;
    if (ch == 0) {
      const int qD = qPi + q1;
      if (qD == 1) {
        if (rng.Uniform() < 2.0 / 3.0) qPiOut = 0, qNOut = 1;
        else qPiOut = 1, qNOut = 0;
      } else if (qD == 0) {
        if (rng.Uniform() < 2.0 / 3.0) qPiOut = 0, qNOut = 0;
        else qPiOut = -1, qNOut = 1;
      }
    }
    const int cPiOut = PionWithCharge(qPiOut), cNOut = NucleonWithCharge(qNOut);
    const double mPiOut = pdg::PionMass(cPiOut), mNOut = pdg::NucleonMass(cNOut);
    if (W <= mPiOut + mNOut + 1e-6) return false;
    const ThreeVector axis = pPi.Boosted(-P.BoostVector()).Vect();
    // The Delta channel's rate already contains the Pauli-blocked width, so its decay is re-drawn
    // until the nucleon is above the Fermi sea; the non-resonant channel is blocked by rejection.
    FourVector p3, p4;
    bool ok = false;
    const int tries = ch == 0 ? 50 : 1;
    for (int it = 0; it < tries && !ok; ++it) {
      double c;
      if (ch == 0) {
        do c = 2 * rng.Uniform() - 1;
        while (4 * rng.Uniform() > 1 + 3 * c * c);  // P33: 1 + 3 cos^2
      } else {
        const double q = hadxs::CMMomentum(W, mPi, pdg::NucleonMass(c1));
        c = SampleExpCos(2 * hadxs::kPiNHighEnergySlope * q * q, rng);
      }
      TwoBody(P, mPiOut, mNOut, AroundAxis(axis, c, 2 * kPi * rng.Uniform()), p3, p4);
      ok = p4.P() >= density.KFermi(cNOut, r);
    }
    if (!ok) {
      if (stats) ++stats->pionPauliBlocked;
      return false;
    }
    if (c1 == pdg::kProton) --rem.z;
    else --rem.n;
    if (stats) {
      if (qPiOut != qPi) ++stats->pionCX;
      else ++stats->pionQE;
    }
    Hadron pi = h;
    pi.pdg = cPiOut;
    pi.dir = p3.Vect().Unit();
    pi.estar = p3.e - mPiOut;
    pi.generation = h.generation + 1;
    active.push_back(pi);
    AddNucleon(nucleus, density, cNOut, h.x, p4, h.ancestor, h.generation + 1, active, rem, stats);
    consumed = true;
    return true;
  }

  // Absorption on two (ch == 2) or three (ch == 3) nucleons. The second nucleon is drawn with the
  // pair isospin weights (np pairs dominate); the third (three-body) from the remnant composition.
  auto pickPartner = [&](int qSoFar, int cFirst, bool pairRule) -> int {
    double w[2];
    for (int s = 0; s < 2; ++s) {
      const int qs = qSoFar + Charge(species[s]);
      const bool allowed = qs >= 0 && qs <= (pairRule ? 2 : 3);
      // availability after the first nucleon(s) are removed
      int have = s == 0 ? rem.z : rem.n;
      if (species[s] == cFirst) have -= 1;
      w[s] = (allowed && have > 0) ? density.RhoSpecies(species[s], r) * have : 0.0;
      if (pairRule && species[s] == cFirst) w[s] *= opt_.likePairAbsorption;
    }
    if (w[0] + w[1] <= 0) return 0;
    return rng.Uniform() * (w[0] + w[1]) < w[0] ? pdg::kProton : pdg::kNeutron;
  };
  const int c2 = pickPartner(qPi + q1, c1, true);
  if (c2 == 0) return false;
  std::vector<int> in = {c1, c2};
  if (ch == 3) {
    // third nucleon: remove c1, c2 from the counts first
    int zA = rem.z - (c1 == pdg::kProton) - (c2 == pdg::kProton), nA = rem.n - (c1 == pdg::kNeutron) - (c2 == pdg::kNeutron);
    const int qSoFar = qPi + q1 + Charge(c2);
    double wp = (qSoFar + 1 <= 3 && zA > 0) ? zA : 0, wn = (qSoFar >= 0 && nA > 0) ? nA : 0;
    if (wp + wn <= 0) return false;
    in.push_back(rng.Uniform() * (wp + wn) < wp ? pdg::kProton : pdg::kNeutron);
  }
  FourVector Ptot = pPi + p2[t];
  int qTot = qPi;
  for (size_t i = 0; i < in.size(); ++i) {
    qTot += Charge(in[i]);
    if (i > 0) Ptot = Ptot + FermiSeaNucleon(in[i], density.KFermi(in[i], r), rng);
  }
  const int nOut = static_cast<int>(in.size());
  if (qTot < 0 || qTot > nOut) return false;
  std::vector<int> outCodes;
  for (int i = 0; i < nOut; ++i) outCodes.push_back(i < qTot ? pdg::kProton : pdg::kNeutron);
  std::vector<double> m(nOut);
  double msum = 0;
  for (int i = 0; i < nOut; ++i) msum += (m[i] = pdg::NucleonMass(outCodes[i]));
  const double M = Ptot.M();
  if (M <= msum + 1e-6) return false;
  std::vector<FourVector> pf(nOut);
  for (int attempt = 0; attempt <= opt_.maxPauliRetries; ++attempt) {
    if (nOut == 2) {
      TwoBody(Ptot, m[0], m[1], IsotropicDirection(rng), pf[0], pf[1]);
    } else {
      // three-body phase space: m23 by rejection on p1* p23*
      const double lo = m[1] + m[2], hi = M - m[0];
      const double wmax = hadxs::CMMomentum(M, m[0], lo) * hadxs::CMMomentum(hi, m[1], m[2]);
      double m23 = lo;
      for (int it = 0; it < 1000; ++it) {
        m23 = rng.Uniform(lo, hi);
        if (rng.Uniform() * wmax <= hadxs::CMMomentum(M, m[0], m23) * hadxs::CMMomentum(m23, m[1], m[2])) break;
      }
      FourVector p23;
      TwoBody(Ptot, m[0], m23, IsotropicDirection(rng), pf[0], p23);
      TwoBody(p23, m[1], m[2], IsotropicDirection(rng), pf[1], pf[2]);
    }
    bool blocked = false;
    for (int i = 0; i < nOut; ++i) blocked |= pf[i].P() < density.KFermi(outCodes[i], r);
    if (!blocked) break;
  }
  for (int c : in) {
    if (c == pdg::kProton) --rem.z;
    else --rem.n;
  }
  if (stats) ++(nOut == 2 ? stats->pionAbs2 : stats->pionAbs3);
  for (int i = 0; i < nOut; ++i)
    AddNucleon(nucleus, density, outCodes[i], h.x, pf[i], h.ancestor, h.generation + 1, active, rem, stats);
  consumed = true;
  return true;
}

}  // namespace nuflow
