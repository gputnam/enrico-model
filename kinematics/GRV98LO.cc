#include "kinematics/GRV98LO.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "core/Config.h"

namespace nuflow {

const GRV98LO& GRV98LO::Instance() {
  static const GRV98LO grv(ResolvePath("data/pdf/grv98lo.dat"));
  return grv;
}

GRV98LO::GRV98LO(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("GRV98LO: cannot open " + path);
  std::string line;
  auto readValues = [&](std::vector<double>& v, size_t n, std::istringstream* first) {
    double d;
    if (first)
      while (v.size() < n && (*first >> d)) v.push_back(d);
    while (v.size() < n && in >> d) v.push_back(d);
    if (v.size() != n) throw std::runtime_error("GRV98LO: truncated table in " + path);
  };
  std::vector<double>* tables[6] = {&uv_, &dv_, &del_, &ud_, &s_, &g_};
  const char* names[6] = {"UV", "DV", "DEL", "UD", "S", "G"};
  while (in >> std::ws && std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ls(line);
    std::string key;
    ls >> key;
    if (key == "NX" || key == "NQ") {
      int n;
      ls >> n;
      if ((key == "NX" && n != kNX) || (key == "NQ" && n != kNQ)) throw std::runtime_error("GRV98LO: bad grid size");
    } else if (key == "LNX") {
      readValues(lnx_, kNX, &ls);
    } else if (key == "LNQ2") {
      readValues(lnq2_, kNQ, &ls);
    } else {
      bool found = false;
      for (int i = 0; i < 6; ++i)
        if (key == names[i]) {
          readValues(*tables[i], kNX * kNQ, nullptr);
          found = true;
        }
      if (!found) throw std::runtime_error("GRV98LO: unknown key '" + key + "' in " + path);
    }
  }
  for (auto* t : tables)
    if (t->size() != static_cast<size_t>(kNX * kNQ)) throw std::runtime_error("GRV98LO: missing table in " + path);
}

namespace {
// Interval and fraction as in the CERNLIB E104 FINT routine: linear inter-/extrapolation between the
// nodes bracketing `a` (the first / last interval outside the grid).
void Locate(const std::vector<double>& nodes, double a, int& i, double& d) {
  const int n = static_cast<int>(nodes.size());
  int j = static_cast<int>(std::lower_bound(nodes.begin(), nodes.end(), a) - nodes.begin());
  if (j >= n) j = n - 1;
  if (j == 0) j = 1;
  i = j - 1;                                            // lower node
  d = (a - nodes[i]) / (nodes[j] - nodes[i]);           // weight of the upper node
}
}  // namespace

double GRV98LO::Interpolate(const std::vector<double>& t, int ix, double dx, int iq, double dq) const {
  auto at = [&](int a, int b) { return t[a + kNX * b]; };
  return (1 - dx) * (1 - dq) * at(ix, iq) + dx * (1 - dq) * at(ix + 1, iq) + (1 - dx) * dq * at(ix, iq + 1) +
         dx * dq * at(ix + 1, iq + 1);
}

PartonDensities GRV98LO::Evaluate(double x, double Q2) const {
  PartonDensities p;
  if (!(x > 0) || x >= 1) return p;
  x = std::max(x, kXMin);
  Q2 = std::clamp(Q2, kQ2Min, kQ2Max);
  int ix, iq;
  double dx, dq;
  Locate(lnx_, std::log(x), ix, dx);
  Locate(lnq2_, std::log(Q2), iq, dq);
  const double x1 = 1 - x, xv = std::sqrt(x), xs = std::pow(x, -0.2);
  const double x13 = x1 * x1 * x1, x14 = x13 * x1, x15 = x14 * x1, x17 = x15 * x1 * x1;
  p.uv = Interpolate(uv_, ix, dx, iq, dq) * x13 * xv;
  p.dv = Interpolate(dv_, ix, dx, iq, dq) * x14 * xv;
  const double del = Interpolate(del_, ix, dx, iq, dq) * x17 * xv;  // x (dbar - ubar)
  const double ud = Interpolate(ud_, ix, dx, iq, dq) * x17 * xs;    // x (ubar + dbar)
  p.ubar = 0.5 * (ud - del);
  p.dbar = 0.5 * (ud + del);
  p.s = Interpolate(s_, ix, dx, iq, dq) * x17 * xs;
  p.g = Interpolate(g_, ix, dx, iq, dq) * x15 * xs;
  return p;
}

}  // namespace nuflow
