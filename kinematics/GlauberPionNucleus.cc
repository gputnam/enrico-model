#include "kinematics/GlauberPionNucleus.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "core/Constants.h"

namespace nuflow {

using namespace constants;

namespace {

constexpr double kHbarC = 0.1973269804;  // GeV fm
constexpr double kGeV2mb = 0.38937937;   // 1 GeV^-2 = 0.3894 mb

// Isoscalar pi N total and elastic cross sections [mb] vs pion lab momentum [GeV]:
// (sigma(pi+ p) + sigma(pi- p))/2 from the PDG data compilation (RPP 2022 files pipp_total, pimp_total,
// pipp_elastic, pimp_elastic; error-weighted means in 0.05 GeV/c bins below 2.5 GeV/c, 0.25 GeV/c to 5 GeV/c,
// logarithmic above), Fermi-averaged in s with a Gaussian of rms 2 p p_F / sqrt(5), p_F = 0.221 GeV (12C).
struct PiNNode {
  double p, tot, el;
};
constexpr PiNNode kPiN[] = {
    {0.45, 33.38, 25.56},   {0.5, 30.74, 22.56},    {0.55, 28.44, 19.62},   {0.6, 27.31, 17.46},
    {0.65, 27.28, 16.13},   {0.7, 27.56, 15.16},    {0.75, 27.87, 14.25},   {0.8, 28.74, 13.64},
    {0.85, 30.56, 13.62},   {0.9, 32.94, 14.19},    {0.95, 35.03, 14.97},   {1, 36.23, 15.56},
    {1.05, 36.52, 15.77},   {1.1, 36.30, 15.65},    {1.15, 36.02, 15.36},   {1.2, 35.94, 15.08},
    {1.25, 36.09, 14.86},   {1.3, 36.38, 14.70},    {1.35, 36.65, 14.54},   {1.4, 36.79, 14.32},
    {1.45, 36.73, 14.03},   {1.5, 36.48, 13.68},    {1.55, 36.08, 13.29},   {1.6, 35.60, 12.88},
    {1.65, 35.08, 12.46},   {1.7, 34.60, 12.05},    {1.75, 34.16, 11.65},   {1.8, 33.80, 11.26},
    {1.85, 33.51, 10.89},   {1.9, 33.28, 10.55},    {1.95, 33.10, 10.23},   {2, 32.97, 9.96},
    {2.05, 32.86, 9.71},    {2.1, 32.78, 9.48},     {2.15, 32.70, 9.28},    {2.2, 32.63, 9.09},
    {2.25, 32.56, 8.91},    {2.3, 32.49, 8.75},     {2.35, 32.40, 8.59},    {2.4, 32.31, 8.44},
    {2.45, 32.21, 8.30},    {2.5, 32.10, 8.18},     {2.75, 31.50, 7.71},    {3, 30.89, 7.43},
    {3.25, 30.38, 7.22},    {3.5, 29.93, 7.01},     {3.75, 29.55, 6.76},    {4, 29.19, 6.50},
    {4.25, 28.85, 6.27},    {4.5, 28.53, 6.07},     {4.75, 28.22, 5.91},    {5, 27.95, 5.79},
    {5.831, 27.25, 5.53},   {6.799, 26.78, 5.30},   {7.929, 26.46, 5.05},   {9.247, 26.06, 4.83},
    {10.78, 25.64, 4.66},   {12.57, 25.29, 4.58},   {14.66, 24.95, 4.40},   {17.1, 24.62, 4.11},
    {19.94, 24.40, 3.82},   {23.25, 24.26, 3.55},   {27.12, 24.11, 3.38},   {31.62, 23.95, 3.36},
    {36.88, 23.82, 3.35},   {43, 23.74, 3.31},      {50.15, 23.70, 3.32},   {58.48, 23.69, 3.34},
    {68.2, 23.70, 3.36},    {79.53, 23.70, 3.34},   {92.74, 23.69, 3.31},   {108.1, 23.70, 3.30},
    {126.1, 23.74, 3.30},   {147.1, 23.81, 3.32},   {171.5, 23.93, 3.33},   {200, 24.01, 3.29},
};
constexpr int kNPiN = sizeof(kPiN) / sizeof(kPiN[0]);

double MandelstamS(double p) {
  const double E = std::sqrt(p * p + kMassPionCharged * kMassPionCharged);
  return kMassProton * kMassProton + kMassPionCharged * kMassPionCharged + 2 * kMassProton * E;
}

// Trapezoid weights for a uniform grid.
double TrapW(int i, int n) { return (i == 0 || i == n - 1) ? 0.5 : 1.0; }

}  // namespace

// Ground-state density of the nucleus.
struct GlauberPionNucleus::Density {
  bool ho = false;
  double a = 0, alpha = 0;  // HO charge density rho ~ (1 + alpha (r/a)^2) exp(-(r/a)^2) [fm]
  double c = 0, d = 0;      // 2pF [fm]
  int A = 0;
  static constexpr double kRN2 = 0.7071 - 0.1161;  // r_p^2 + r_n^2 [fm^2] (CODATA 2018 / PDG)

  static Density For(int A) {
    Density D;
    D.A = A;
    // de Vries, de Jager, de Vries, At. Data Nucl. Data Tables 36, 495 (1987).
    if (A == 12) { D.ho = true; D.a = 1.687; D.alpha = 1.067; return D; }
    if (A == 16) { D.ho = true; D.a = 1.833; D.alpha = 1.544; return D; }
    if (A == 27) { D.c = 3.07; D.d = 0.519; return D; }
    if (A == 28) { D.c = 3.14; D.d = 0.537; return D; }
    if (A == 40) { D.c = 3.53; D.d = 0.542; return D; }
    if (A == 56) { D.c = 4.106; D.d = 0.519; return D; }
    // generic 2pF: rms radius from the global charge-radius systematics, d = 0.537 fm
    const double a13 = std::cbrt(double(A));
    const double rrms = 0.891 * a13 * (1 + 1.565 / (a13 * a13) - 1.043 / (a13 * a13 * a13 * a13));
    D.d = 0.537;
    D.c = std::sqrt(std::max(0.25, (5 * rrms * rrms - 7 * kPi * kPi * D.d * D.d) / 3));
    return D;
  }
  // HO point-nucleon oscillator width used for the c.m. (fixed-centre) correction [fm].
  double bHO2() const { return a * a - 2.0 * kRN2 / 3.0; }
};

struct GlauberPionNucleus::Cache {
  std::mutex m;
  std::array<std::atomic<const Table*>, 301> tables{};
  std::vector<std::unique_ptr<Table>> owned;
};

GlauberPionNucleus::GlauberPionNucleus() : cache_(std::make_shared<Cache>()) {
  for (auto& t : cache_->tables) t.store(nullptr);
}

double GlauberPionNucleus::ReggeRho(double p) {
  // PDG Regge + ln^2 s fit of the isoscalar Im f (Z, B ln^2(s/sM), Y1 (s1/s)^eta1) and its derivative-
  // analyticity real part: Re = -cot(pi alpha/2) * Y1 s^-eta1 (C = +1 Regge, alpha = 1 - eta1) + pi B ln(s/sM).
  constexpr double Z = 20.86, Y1 = 19.24, eta1 = 0.4473, M = 2.1206;
  const double B = kPi * kGeV2mb / (M * M);
  const double sM = std::pow(kMassPionCharged + kMassProton + M, 2);
  const double s = MandelstamS(p);
  const double l = std::log(s / sM);
  const double reg = Y1 * std::pow(s, -eta1);
  const double im = Z + B * l * l + reg;
  const double re = -reg / std::tan(kPi * (1 - eta1) / 2) + kPi * B * l;
  return re / im;
}

GlauberPionNucleus::PiN GlauberPionNucleus::PiNAmplitude(double p) {
  PiN r;
  const double pc = std::clamp(p, kPiN[0].p, kPiN[kNPiN - 1].p);
  int i = 0;
  while (i < kNPiN - 2 && pc > kPiN[i + 1].p) ++i;
  const double f = std::log(pc / kPiN[i].p) / std::log(kPiN[i + 1].p / kPiN[i].p);
  r.sigmaTot = kPiN[i].tot + f * (kPiN[i + 1].tot - kPiN[i].tot);
  r.sigmaEl = kPiN[i].el + f * (kPiN[i + 1].el - kPiN[i].el);
  r.rho = ReggeRho(p);
  const double st = r.sigmaTot / kGeV2mb, se = r.sigmaEl / kGeV2mb;  // GeV^-2
  r.beta = st * st * (1 + r.rho * r.rho) / (16 * kPi * se);
  return r;
}

std::vector<double> GlauberPionNucleus::Thickness(double T, int A, PiN& pin) const {
  const Density D = Density::For(A);
  const double p = std::sqrt(T * T + 2 * T * kMassPionCharged);
  pin = PiNAmplitude(p);
  const double betaFm = pin.beta * kHbarC * kHbarC;  // fm^2
  const double db = kBMax / (kNb - 1);
  std::vector<double> Tb(kNb);
  if (D.ho) {
    // S(q) e^{-beta q^2/2} = (1 - c q^2) e^{-g q^2}: analytic Hankel transform.
    const double c = D.alpha * D.a * D.a / (2 * (2 + 3 * D.alpha));
    const double g = D.a * D.a / 4 - Density::kRN2 / 6 + D.bHO2() / (4.0 * A) + betaFm / 2;
    for (int i = 0; i < kNb; ++i) {
      const double b = i * db, x = b * b / (4 * g), e = std::exp(-x);
      Tb[i] = (e / (2 * g) - c * e * (1 - x) / (2 * g * g)) / (2 * kPi);
    }
  } else {
    // 2pF (point density): numerical 3D -> 2D: T(b) = Int dz rho, then Gaussian folding with the piN range
    // done in q space via the numerical form factor.
    constexpr int nr = 600, nq = 400;
    const double rmax = D.c + 15 * D.d, dr = rmax / (nr - 1), qmax = 8.0, dq = qmax / (nq - 1);
    std::vector<double> rho(nr), S(nq);
    double norm = 0;
    for (int i = 0; i < nr; ++i) {
      const double r = i * dr;
      rho[i] = 1.0 / (1 + std::exp((r - D.c) / D.d));
      norm += TrapW(i, nr) * 4 * kPi * r * r * rho[i] * dr;
    }
    for (int j = 0; j < nq; ++j) {
      const double q = j * dq;
      double s = 0;
      for (int i = 0; i < nr; ++i) {
        const double r = i * dr, x = q * r;
        s += TrapW(i, nr) * 4 * kPi * r * r * rho[i] * (x > 1e-8 ? std::sin(x) / x : 1.0) * dr;
      }
      S[j] = s / norm * std::exp(-betaFm * q * q / 2);
    }
    for (int i = 0; i < kNb; ++i) {
      const double b = i * db;
      double s = 0;
      for (int j = 0; j < nq; ++j) s += TrapW(j, nq) * std::cyl_bessel_j(0.0, j * dq * b) * S[j] * j * dq * dq;
      Tb[i] = s / (2 * kPi);
    }
  }
  return Tb;
}

std::vector<std::complex<double>> GlauberPionNucleus::Profile(double T, int A, PiN* used) const {
  PiN pin;
  const std::vector<double> Tb = Thickness(T, A, pin);
  if (used) *used = pin;
  const double sigFm = pin.sigmaTot / 10.0;  // fm^2
  std::vector<std::complex<double>> G(kNb);
  const std::complex<double> x0 = std::complex<double>(1.0, -pin.rho) * (sigFm / 2);
  for (int i = 0; i < kNb; ++i) G[i] = 1.0 - std::pow(1.0 - x0 * Tb[i], A);
  return G;
}

GlauberPionNucleus::Observables GlauberPionNucleus::Compute(double T, int A) const {
  PiN pin;
  const auto G = Profile(T, A, &pin);
  const auto Tb = Thickness(T, A, pin);
  const double db = kBMax / (kNb - 1);
  const double sIn = (pin.sigmaTot - pin.sigmaEl) / 10.0;
  Observables o;
  for (int i = 0; i < kNb; ++i) {
    const double b = i * db, w = TrapW(i, kNb) * 2 * kPi * b * db * 10.0;  // fm^2 -> mb
    o.sigmaTot += 2 * G[i].real() * w;
    o.sigmaEl += std::norm(G[i]) * w;
    // production ("absorption"): Int d2b {1 - [1 - sigma_in T(b)]^A}, sigma_in = sigma_tot - sigma_el (piN)
    o.sigmaProd += (1 - std::pow(std::max(0.0, 1 - sIn * Tb[i]), A)) * w;
  }
  o.sigmaInel = o.sigmaTot - o.sigmaEl;
  return o;
}

double GlauberPionNucleus::DSigmaDtDirect(double T, double absT, int A) const {
  const auto G = Profile(T, A);
  const Density D = Density::For(A);
  const double db = kBMax / (kNb - 1);
  const double q = std::sqrt(std::max(0.0, absT)) / kHbarC;  // fm^-1
  std::complex<double> I = 0;
  for (int i = 0; i < kNb; ++i) {
    const double b = i * db;
    I += TrapW(i, kNb) * b * std::cyl_bessel_j(0.0, q * b) * G[i] * db;
  }
  if (D.ho) I *= std::exp(q * q * D.bHO2() / (4.0 * A));
  return kPi * std::norm(I) / (kHbarC * kHbarC) * 10.0;  // fm^4 -> mb/GeV^2
}

const GlauberPionNucleus::Table& GlauberPionNucleus::GetTable(int A) const {
  if (A < 1 || A > 300) throw std::runtime_error("GlauberPionNucleus: A out of range");
  if (const Table* t = cache_->tables[A].load(std::memory_order_acquire)) return *t;
  std::lock_guard<std::mutex> lock(cache_->m);
  if (const Table* t = cache_->tables[A].load(std::memory_order_acquire)) return *t;
  auto tab = std::make_unique<Table>();
  const Density D = Density::For(A);
  const double db = kBMax / (kNb - 1), dq = kQMax / (kNq - 1);
  // J0(q b) b db weights
  std::vector<double> J(size_t(kNq) * kNb);
  for (int j = 0; j < kNq; ++j)
    for (int i = 0; i < kNb; ++i) {
      const double q = j * dq / kHbarC, b = i * db;
      J[size_t(j) * kNb + i] = TrapW(i, kNb) * b * db * std::cyl_bessel_j(0.0, q * b);
    }
  tab->lnT.resize(kNT);
  tab->amp.assign(kNT, std::vector<std::complex<double>>(kNq));
  for (int k = 0; k < kNT; ++k) {
    const double lnT = std::log(kTMin) + k * (std::log(kTMax) - std::log(kTMin)) / (kNT - 1);
    tab->lnT[k] = lnT;
    const auto G = Profile(std::exp(lnT), A);
    for (int j = 0; j < kNq; ++j) {
      std::complex<double> I = 0;
      const double* Jj = &J[size_t(j) * kNb];
      for (int i = 0; i < kNb; ++i) I += Jj[i] * G[i];
      const double q = j * dq / kHbarC;
      if (D.ho) I *= std::exp(q * q * D.bHO2() / (4.0 * A));
      tab->amp[k][j] = I;
    }
  }
  const Table* raw = tab.get();
  cache_->owned.push_back(std::move(tab));
  cache_->tables[A].store(raw, std::memory_order_release);
  return *raw;
}

double GlauberPionNucleus::DSigmaDt(double T, double absT, int A) const {
  const Table& t = GetTable(A);
  const double q = std::sqrt(std::max(0.0, absT));
  if (q >= kQMax) return 0;
  const double lnT = std::log(std::clamp(T, kTMin, kTMax));
  const double fk = (lnT - t.lnT[0]) / (t.lnT[1] - t.lnT[0]);
  const int k = std::clamp(int(fk), 0, kNT - 2);
  const double wk = std::clamp(fk - k, 0.0, 1.0);
  const double fq = q / (kQMax / (kNq - 1));
  const int j = std::clamp(int(fq), 0, kNq - 2);
  const double wq = fq - j;
  const std::complex<double> I = (1 - wk) * ((1 - wq) * t.amp[k][j] + wq * t.amp[k][j + 1]) +
                                 wk * ((1 - wq) * t.amp[k + 1][j] + wq * t.amp[k + 1][j + 1]);
  return kPi * std::norm(I) / (kHbarC * kHbarC) * 10.0;  // mb/GeV^2
}

}  // namespace nuflow
