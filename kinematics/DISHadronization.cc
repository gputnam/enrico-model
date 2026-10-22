#include "kinematics/DISHadronization.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <mutex>
#include <numeric>

#include "core/Constants.h"
#include "core/PDG.h"

namespace nuflow {

using namespace constants;

namespace {

double Normal(Random& rng) {
  const double u1 = rng.UniformOpen(), u2 = rng.Uniform();
  return std::sqrt(-2 * std::log(u1)) * std::cos(2 * kPi * u2);
}

// Marsaglia-Tsang gamma variate with shape k and unit scale.
double Gamma(double k, Random& rng) {
  if (k < 1) return Gamma(k + 1, rng) * std::pow(rng.UniformOpen(), 1 / k);
  const double d = k - 1.0 / 3, c = 1 / std::sqrt(9 * d);
  for (;;) {
    double z, v;
    do {
      z = Normal(rng);
      v = 1 + c * z;
    } while (v <= 0);
    v = v * v * v;
    const double u = rng.UniformOpen();
    if (std::log(u) < 0.5 * z * z + d - d * v + d * std::log(v)) return d * v;
  }
}

int Poisson(double mean, Random& rng) {
  if (mean <= 0) return 0;
  const double L = std::exp(-mean);
  int k = 0;
  double p = rng.Uniform();
  while (p > L) {
    ++k;
    p *= rng.Uniform();
  }
  return k;
}

double Kallen(double a, double b, double c) { return a * a + b * b + c * c - 2 * (a * b + a * c + b * c); }
double TwoBodyP(double M, double m1, double m2) {
  const double l = Kallen(M * M, m1 * m1, m2 * m2);
  return l > 0 ? std::sqrt(l) / (2 * M) : 0;
}

ThreeVector IsotropicDirection(Random& rng) {
  const double c = 2 * rng.Uniform() - 1;
  const double s = std::sqrt(std::max(0.0, 1 - c * c));
  const double phi = 2 * kPi * rng.Uniform();
  return {s * std::cos(phi), s * std::sin(phi), c};
}

double Mass(int code) {
  if (code == pdg::kProton || code == pdg::kNeutron) return pdg::NucleonMass(code);
  return pdg::PionMass(code);
}

int NumPions(const std::vector<int>& pdgs) {
  int n = 0;
  for (int c : pdgs)
    if (c != pdg::kProton && c != pdg::kNeutron) ++n;
  return n;
}

}  // namespace

DISHadronization::DISHadronization(Options o) : opt_(o) { BuildSinglePionTable(); }

void DISHadronization::BuildSinglePionTable() {
  // P_1pi(W) by Monte Carlo with a fixed seed: 500 nodes in 1.07 < W < 6.06 GeV, 10000 draws per node
  // (binomial error <= 0.005). Only the dispersion slope enters the species sampling, so tables are
  // shared between instances with the same slope.
  static std::mutex mutex;
  static std::map<double, std::shared_ptr<const P1Table>> cache;
  std::lock_guard<std::mutex> lock(mutex);
  if (auto it = cache.find(opt_.dispersionSlope); it != cache.end()) {
    p1_ = it->second;
    return;
  }
  auto t = std::make_shared<P1Table>();
  t->wLo = 1.07;
  t->dW = 0.01;
  t->n = 500;
  constexpr int kDraws = 10000;
  const int nus[2] = {14, -14};
  const int hits[2] = {pdg::kProton, pdg::kNeutron};
  std::vector<int> pdgs;
  for (int a = 0; a < 2; ++a)
    for (int b = 0; b < 2; ++b) {
      auto& v = t->p[2 * a + b];
      v.resize(t->n);
      Random rng(0x5151A7ULL, 2 * a + b);
      for (int i = 0; i < t->n; ++i) {
        const double W = t->wLo + i * t->dW;
        if (W <= kMassNeutron + kMassPionCharged) {
          v[i] = 1;
          continue;
        }
        int single = 0;
        for (int k = 0; k < kDraws; ++k) {
          SampleSpecies(nus[a], hits[b], W, rng, pdgs);
          if (NumPions(pdgs) == 1) ++single;
        }
        v[i] = static_cast<double>(single) / kDraws;
      }
    }
  p1_ = t;
  cache[opt_.dispersionSlope] = t;
}

double DISHadronization::SinglePionFraction(int nuPdg, int hitPdg, double W) const {
  const auto& v = p1_->p[2 * (nuPdg < 0 ? 1 : 0) + (hitPdg == pdg::kNeutron ? 1 : 0)];
  const double u = (W - p1_->wLo) / p1_->dW;
  if (u <= 0) return v.front();
  if (u >= p1_->n - 1) return v.back();
  const int i = static_cast<int>(u);
  const double f = u - i;
  return (1 - f) * v[i] + f * v[i + 1];
}

bool DISHadronization::SampleSpecies(int nuPdg, int hitPdg, double W, PionClass cls, Random& rng,
                                     std::vector<int>& pdgs) const {
  if (cls == PionClass::kAny) {
    SampleSpecies(nuPdg, hitPdg, W, rng, pdgs);
    return true;
  }
  const double mMin2pi = kMassProton + 2 * kMassPionNeutral;
  if (cls == PionClass::kMulti && W <= mMin2pi + 1e-3) return false;
  for (int attempt = 0; attempt < 20000; ++attempt) {
    SampleSpecies(nuPdg, hitPdg, W, rng, pdgs);
    const int n = NumPions(pdgs);
    if ((cls == PionClass::kSingle) == (n == 1)) return true;
  }
  return false;
}

double DISHadronization::MeanChargedMultiplicity(int nuPdg, int hitPdg, double W) {
  // <n_ch> = a + b ln W^2, deuterium bubble chambers (see header).
  double a, b;
  const bool proton = hitPdg == pdg::kProton;
  if (nuPdg > 0) {
    a = proton ? 0.50 : -0.20;  // FNAL E545 (Zieminska 1983), nu p / nu n
    b = 1.42;
  } else {
    a = proton ? 0.02 : 0.80;  // CERN WA25 (Barlag 1982), nubar p / nubar n
    b = proton ? 1.28 : 0.95;
  }
  return a + b * std::log(W * W);
}

void DISHadronization::SampleSpecies(int nuPdg, int hitPdg, double W, Random& rng, std::vector<int>& pdgs) const {
  const bool nu = nuPdg > 0;
  const bool proton = hitPdg == pdg::kProton;
  const int QX = (proton ? 1 : 0) + (nu ? 1 : -1);  // hadronic charge
  const double nch = MeanChargedMultiplicity(nuPdg, hitPdg, W);
  const double muPair = std::max(0.0, 0.5 * (nch - std::abs(QX)));
  const double sigPair = 0.5 * opt_.dispersionSlope * std::max(nch, 0.0);
  // Variance of the continuous variate, leaving room for the stochastic rounding (~1/6 on average).
  const double var = std::max(sigPair * sigPair - 1.0 / 6, 0.09 * sigPair * sigPair);

  for (int attempt = 0; attempt < 1000; ++attempt) {
    pdgs.clear();
    // Leading baryon: spectator diquark (uu for nu p, ud for nu n / nubar p, dd for nubar n) plus
    // a u or d quark with equal probability (uuu -> Delta++ -> p pi+, ddd -> Delta- -> n pi-).
    int baryon;
    if (nu) baryon = proton ? pdg::kProton : (rng.Uniform() < 0.5 ? pdg::kProton : pdg::kNeutron);
    else baryon = proton ? (rng.Uniform() < 0.5 ? pdg::kProton : pdg::kNeutron) : pdg::kNeutron;
    const int dp = baryon == pdg::kProton ? 1 : 0;
    const int QM = QX - dp;  // meson charge

    // Number of +- pairs beyond the minimum charged multiplicity |Q_X|.
    int nPair = 0;
    if (muPair > 0) {
      const double k = muPair * muPair / var, theta = var / muPair;
      const double t = Gamma(k, rng) * theta;
      nPair = static_cast<int>(std::floor(t));
      if (rng.Uniform() < t - nPair) ++nPair;
    }
    int nch_ = std::abs(QX) + 2 * nPair;
    while (nch_ < dp + std::abs(QM)) nch_ += 2;  // charge conservation with a charged baryon
    const int nPip = (nch_ - dp + QM) / 2, nPim = (nch_ - dp - QM) / 2;
    int nPi0 = Poisson(0.5 * (nPip + nPim), rng);
    if (nPip + nPim + nPi0 == 0) nPi0 = 1;  // W > M + m_pi: at least one meson

    pdgs.push_back(baryon);
    for (int i = 0; i < nPip; ++i) pdgs.push_back(pdg::kPiPlus);
    for (int i = 0; i < nPim; ++i) pdgs.push_back(pdg::kPiMinus);
    for (int i = 0; i < nPi0; ++i) pdgs.push_back(pdg::kPiZero);
    double msum = 0;
    for (int c : pdgs) msum += Mass(c);
    if (msum < W - 1e-3) return;
  }
  // Fallback (never reached for W above the two-body threshold): minimal final state.
  pdgs.clear();
  if (nu) {
    pdgs = proton ? std::vector<int>{pdg::kProton, pdg::kPiPlus} : std::vector<int>{pdg::kProton, pdg::kPiZero};
  } else {
    pdgs = proton ? std::vector<int>{pdg::kNeutron, pdg::kPiZero} : std::vector<int>{pdg::kNeutron, pdg::kPiMinus};
  }
}

bool DISHadronization::PhaseSpace(double W, const std::vector<double>& m, Random& rng, std::vector<FourVector>& p) {
  // Raubold-Lynch / GENBOD (F. James, CERN 68-15): ordered random invariant masses, weight = product
  // of two-body momenta, unweighted by rejection against an upper bound.
  const int n = static_cast<int>(m.size());
  p.assign(n, FourVector());
  const double msum = std::accumulate(m.begin(), m.end(), 0.0);
  const double T = W - msum;
  if (T <= 0) return false;
  if (n == 1) {
    p[0] = FourVector(W, 0, 0, 0);
    return true;
  }
  // Upper bound of the weight.
  double wmax = 1;
  {
    double emmax = T + m[0], emmin = 0;
    for (int i = 1; i < n; ++i) {
      emmin += m[i - 1];
      emmax += m[i];
      wmax *= TwoBodyP(emmax, emmin, m[i]);
    }
  }
  std::vector<double> r(n), M(n), pd(n);
  for (int attempt = 0; attempt < 100000; ++attempt) {
    r[0] = 0;
    r[n - 1] = 1;
    for (int i = 1; i < n - 1; ++i) r[i] = rng.Uniform();
    std::sort(r.begin() + 1, r.end() - 1);
    double sum = 0, w = 1;
    for (int i = 0; i < n; ++i) {
      sum += m[i];
      M[i] = r[i] * T + sum;
    }
    for (int i = 0; i + 1 < n; ++i) {
      pd[i] = TwoBodyP(M[i + 1], M[i], m[i + 1]);
      w *= pd[i];
    }
    if (rng.Uniform() * wmax > w) continue;
    // Build momenta: successive two-body decays M[i+1] -> M[i] + m[i+1], each along an independent
    // isotropic direction; the subsystem 0..i is boosted into the rest frame of M[i+1].
    const ThreeVector d0 = IsotropicDirection(rng);
    p[0] = FourVector(std::sqrt(pd[0] * pd[0] + m[0] * m[0]), d0 * pd[0]);
    p[1] = FourVector(std::sqrt(pd[0] * pd[0] + m[1] * m[1]), d0 * (-pd[0]));
    for (int i = 1; i + 1 < n; ++i) {
      const ThreeVector d = IsotropicDirection(rng);
      const double e = std::sqrt(pd[i] * pd[i] + M[i] * M[i]);
      const ThreeVector beta = d * (-pd[i] / e);  // subsystem 0..i recoils along -d
      for (int j = 0; j <= i; ++j) p[j] = p[j].Boosted(beta);
      p[i + 1] = FourVector(std::sqrt(pd[i] * pd[i] + m[i + 1] * m[i + 1]), d * pd[i]);
    }
    return true;
  }
  return false;
}

bool DISHadronization::LongitudinalPhaseSpace(double W, const std::vector<double>& m, int backward, Random& rng,
                                              std::vector<FourVector>& p) const {
  const int n = static_cast<int>(m.size());
  p.assign(n, FourVector());
  std::vector<double> px(n), py(n), mt(n), y(n);
  const double sig = std::sqrt(0.5 * opt_.meanPt2);
  for (int attempt = 0; attempt < 1000; ++attempt) {
    // Transverse momenta: Gaussian, then shifted so that sum pT = 0. After many failures the pT
    // scale is reduced (only relevant just above threshold).
    const double scale = attempt < 200 ? 1.0 : std::pow(0.98, attempt - 199);
    double sx = 0, sy = 0;
    for (int i = 0; i < n; ++i) {
      px[i] = sig * scale * Normal(rng);
      py[i] = sig * scale * Normal(rng);
      sx += px[i];
      sy += py[i];
    }
    double mtsum = 0;
    for (int i = 0; i < n; ++i) {
      px[i] -= sx / n;
      py[i] -= sy / n;
      mt[i] = std::sqrt(m[i] * m[i] + px[i] * px[i] + py[i] * py[i]);
      mtsum += mt[i];
    }
    if (mtsum >= W) continue;
    // Rapidities: uniform plateau; the baryon gets the smallest (target fragmentation side).
    std::vector<double> u(n);
    for (auto& v : u) v = rng.Uniform();
    std::sort(u.begin(), u.end());
    // mesons take the remaining values in random order
    std::vector<int> order;
    for (int i = 0; i < n; ++i)
      if (i != backward) order.push_back(i);
    for (int i = static_cast<int>(order.size()) - 1; i > 0; --i)
      std::swap(order[i], order[static_cast<int>(rng.Uniform() * (i + 1))]);
    if (backward >= 0) {
      y[backward] = u[0];
      for (size_t j = 0; j < order.size(); ++j) y[order[j]] = u[j + 1];
    } else {
      for (size_t j = 0; j < order.size(); ++j) y[order[j]] = u[j];
    }
    if (u.back() - u.front() < 1e-3) continue;
    // Stretch the rapidities by c and centre them so that the invariant mass is W and pz = 0.
    // M^2(c) = sum_ij mT_i mT_j cosh(c (y_i - y_j)) is free of cancellations and increases with c.
    auto mass2 = [&](double c) {
      double m2 = 0;
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) m2 += mt[i] * mt[j] * std::cosh(c * (y[i] - y[j]));
      return m2;
    };
    const double W2 = W * W;
    double lo = 0, hi = 1;
    while (mass2(hi) < W2 && hi < 1e4) hi *= 2;
    if (mass2(hi) < W2) continue;
    for (int it = 0; it < 200 && hi - lo > 1e-15 * hi; ++it) {
      const double mid = 0.5 * (lo + hi);
      (mass2(mid) < W2 ? lo : hi) = mid;
    }
    const double c = 0.5 * (lo + hi);
    // Rapidity of the system, Y = 0.5 ln[(E + P)/(E - P)], via log-sum-exp.
    double ymax = -1e300, ymin = 1e300;
    for (int i = 0; i < n; ++i) {
      ymax = std::max(ymax, c * y[i]);
      ymin = std::min(ymin, c * y[i]);
    }
    double sp = 0, sm = 0;
    for (int i = 0; i < n; ++i) {
      sp += mt[i] * std::exp(c * y[i] - ymax);
      sm += mt[i] * std::exp(-(c * y[i] - ymin));
    }
    const double Y = 0.5 * ((ymax + std::log(sp)) - (-ymin + std::log(sm)));
    for (int i = 0; i < n; ++i) {
      const double yy = c * y[i] - Y;
      p[i] = FourVector(mt[i] * std::cosh(yy), px[i], py[i], mt[i] * std::sinh(yy));
    }
    return true;
  }
  return false;
}

void DISHadronization::Hadronize(int nuPdg, int hitPdg, const FourVector& X, const ThreeVector& axis, Random& rng,
                                 std::vector<Hadron>& out, PionClass cls) const {
  HadronizeImpl(nuPdg, hitPdg, X, axis, rng, out, cls, nullptr);
}

void DISHadronization::HadronizeSpecies(const std::vector<int>& species, const FourVector& X, const ThreeVector& axis,
                                        Random& rng, std::vector<Hadron>& out) const {
  HadronizeImpl(0, 0, X, axis, rng, out, PionClass::kAny, &species);
}

void DISHadronization::HadronizeImpl(int nuPdg, int hitPdg, const FourVector& X, const ThreeVector& axis, Random& rng,
                                     std::vector<Hadron>& out, PionClass cls, const std::vector<int>* fixed) const {
  out.clear();
  const double W = X.M();
  std::vector<int> pdgs;
  std::vector<double> masses;
  std::vector<FourVector> p;
  for (int attempt = 0; attempt < 100; ++attempt) {
    if (fixed) pdgs = *fixed;
    else if (!SampleSpecies(nuPdg, hitPdg, W, cls, rng, pdgs)) return;
    masses.clear();
    for (int c : pdgs) masses.push_back(Mass(c));
    const bool ok = (W < opt_.wLPS || pdgs.size() == 2) ? PhaseSpace(W, masses, rng, p)
                                                         : LongitudinalPhaseSpace(W, masses, 0, rng, p);
    if (ok) break;
    p.clear();
  }
  if (p.size() != pdgs.size()) return;  // caller treats an empty list as a failure
  // Rest frame of X with z along the momentum transfer (direction of `axis` boosted into it).
  const ThreeVector beta = X.BoostVector();
  const FourVector ax = FourVector(axis.Mag(), axis).Boosted(-beta);
  const ThreeVector ez = ax.Vect().Unit();
  ThreeVector ex = (std::abs(ez.x) < 0.9 ? ThreeVector{1, 0, 0} : ThreeVector{0, 1, 0}).Cross(ez).Unit();
  ThreeVector ey = ez.Cross(ex);
  // random azimuth about the axis
  const double phi = 2 * kPi * rng.Uniform();
  const ThreeVector e1 = ex * std::cos(phi) + ey * std::sin(phi), e2 = ey * std::cos(phi) - ex * std::sin(phi);
  for (size_t i = 0; i < pdgs.size(); ++i) {
    const ThreeVector v = p[i].Vect();
    const FourVector rest(p[i].e, e1 * v.x + e2 * v.y + ez * v.z);
    out.push_back({pdgs[i], rest.Boosted(beta)});
  }
}

}  // namespace nuflow
