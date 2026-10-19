#include "orchestration/EnergySampler.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace nuflow {

EnergySampler::EnergySampler(const FluxSource& flux, const ChannelSelector& selector) : selector_(selector) {
  if (flux.IsMonoEnergetic()) {
    mono_ = true;
    monoEnergy_ = flux.EMin();
    fluxAveragedSigma_ = selector.TotalSigma(monoEnergy_);
    if (fluxAveragedSigma_ <= 0) throw std::runtime_error("EnergySampler: zero cross section at beam energy");
    return;
  }
  const auto& hist = dynamic_cast<const HistogramFlux&>(flux);
  edges_ = hist.Edges();
  const auto& dens = hist.Density();
  // All table nodes (union over processes) to locate the maxima of the piecewise-linear sigma(E).
  std::vector<double> nodes;
  for (size_t i = 0; i < selector.Size(); ++i) {
    const auto& e = selector.Get(i).table.energy;
    nodes.insert(nodes.end(), e.begin(), e.end());
  }
  std::sort(nodes.begin(), nodes.end());

  const size_t nb = dens.size();
  cdf_.assign(nb, 0);
  sigmaMax_.assign(nb, 0);
  double total = 0, fluxIntegral = 0;
  for (size_t i = 0; i < nb; ++i) {
    const double lo = edges_[i], hi = edges_[i + 1];
    // Integrate sigma over the bin exactly for a piecewise-linear interpolant: trapezoid between
    // consecutive breakpoints (bin edges + interior nodes).
    std::vector<double> pts{lo};
    for (auto it = std::upper_bound(nodes.begin(), nodes.end(), lo); it != nodes.end() && *it < hi; ++it)
      pts.push_back(*it);
    pts.push_back(hi);
    double integral = 0, smax = 0;
    for (size_t j = 0; j < pts.size(); ++j) {
      const double s = selector.TotalSigma(pts[j]);
      smax = std::max(smax, s);
      if (j > 0) integral += 0.5 * (s + selector.TotalSigma(pts[j - 1])) * (pts[j] - pts[j - 1]);
    }
    // Guard against the step at the first table node (sigma = 0 below it).
    smax = std::max(smax, selector.TotalSigma(std::nextafter(hi, lo)));
    sigmaMax_[i] = smax;
    total += dens[i] * integral;
    fluxIntegral += dens[i] * (hi - lo);
    cdf_[i] = total;
  }
  if (total <= 0) throw std::runtime_error("EnergySampler: flux x cross section integrates to zero");
  for (auto& c : cdf_) c /= total;
  fluxAveragedSigma_ = total / fluxIntegral;
}

double EnergySampler::Sample(Random& rng) const {
  if (mono_) return monoEnergy_;
  const size_t i = std::min<size_t>(std::upper_bound(cdf_.begin(), cdf_.end(), rng.Uniform()) - cdf_.begin(),
                                    cdf_.size() - 1);
  const double lo = edges_[i], hi = edges_[i + 1];
  // Within the bin the flux is constant: sample E ~ sigma(E) by rejection.
  for (;;) {
    const double e = lo + (hi - lo) * rng.Uniform();
    if (rng.Uniform() * sigmaMax_[i] <= selector_.TotalSigma(e)) return e;
  }
}

}  // namespace nuflow
