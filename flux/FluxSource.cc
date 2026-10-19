#include "flux/FluxSource.h"

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace nuflow {

std::string MonoEnergeticFlux::Description() const {
  std::ostringstream s;
  s << "mono-energetic " << energy_ << " GeV";
  return s.str();
}

HistogramFlux::HistogramFlux(std::vector<double> edges, std::vector<double> density, std::string description,
                             std::string units)
    : edges_(std::move(edges)), density_(std::move(density)), description_(std::move(description)),
      units_(std::move(units)) {
  if (edges_.size() != density_.size() + 1 || density_.empty())
    throw std::runtime_error("HistogramFlux: need N+1 edges for N bins");
  for (size_t i = 0; i < density_.size(); ++i) {
    if (!(edges_[i + 1] > edges_[i])) throw std::runtime_error("HistogramFlux: bin edges must increase");
    if (density_[i] < 0) throw std::runtime_error("HistogramFlux: negative flux in bin " + std::to_string(i));
  }
}

double HistogramFlux::Integral() const {
  double sum = 0;
  for (size_t i = 0; i < density_.size(); ++i) sum += density_[i] * (edges_[i + 1] - edges_[i]);
  return sum;
}

void HistogramFlux::Clip(double emin, double emax) {
  std::vector<double> e, d;
  for (size_t i = 0; i < density_.size(); ++i) {
    const double lo = std::max(edges_[i], emin), hi = std::min(edges_[i + 1], emax);
    if (hi <= lo) continue;
    if (e.empty()) e.push_back(lo);
    else if (e.back() != lo) throw std::runtime_error("HistogramFlux::Clip: non-contiguous bins");
    e.push_back(hi);
    d.push_back(density_[i]);
  }
  if (d.empty()) throw std::runtime_error("HistogramFlux::Clip: empty energy range");
  // Drop leading/trailing empty bins so the generator's energy range is tight.
  size_t first = 0, last = d.size();
  while (first < last && d[first] == 0) ++first;
  while (last > first && d[last - 1] == 0) --last;
  if (first == last) throw std::runtime_error("HistogramFlux: flux is zero everywhere");
  edges_.assign(e.begin() + first, e.begin() + last + 1);
  density_.assign(d.begin() + first, d.begin() + last);
}

}  // namespace nuflow
