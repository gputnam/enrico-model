#include "channel/ChannelSelector.h"

#include <algorithm>
#include <stdexcept>

namespace nuflow {

size_t XSecTable::Bracket(double e) const {
  // index i such that energy[i] <= e < energy[i+1] (callers handle the ends)
  return static_cast<size_t>(std::upper_bound(energy.begin(), energy.end(), e) - energy.begin()) - 1;
}

double XSecTable::Sigma(double e) const {
  if (energy.size() == 1) return sigma[0];
  if (e <= energy.front()) return e == energy.front() ? sigma.front() : 0.0;
  if (e >= energy.back()) return sigma.back();
  const size_t i = Bracket(e);
  const double t = (e - energy[i]) / (energy[i + 1] - energy[i]);
  return sigma[i] + t * (sigma[i + 1] - sigma[i]);
}

double XSecTable::WMax(double e) const {
  if (energy.size() == 1 || e <= energy.front()) return wMax.front();
  if (e >= energy.back()) return wMax.back();
  const size_t i = Bracket(e);
  return std::max(wMax[i], wMax[i + 1]);
}

double ChannelSelector::TotalSigma(double e) const {
  double s = 0;
  for (const auto& p : processes_) s += p.count * p.table.Sigma(e);
  return s;
}

size_t ChannelSelector::Select(double e, Random& rng) const {
  // Few processes: a linear cumulative scan is faster than building alias tables per energy.
  double cum[64];
  const size_t n = processes_.size();
  if (n > 64) throw std::logic_error("ChannelSelector: too many processes");
  double total = 0;
  for (size_t i = 0; i < n; ++i) cum[i] = (total += processes_[i].count * processes_[i].table.Sigma(e));
  if (total <= 0) throw std::runtime_error("ChannelSelector: zero cross section at this energy");
  const double u = rng.Uniform() * total;
  for (size_t i = 0; i < n; ++i)
    if (u < cum[i]) return i;
  return n - 1;
}

}  // namespace nuflow
