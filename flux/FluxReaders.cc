// Readers for flux files. All produce a HistogramFlux with densities per GeV.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

#include <TFile.h>
#include <TH1.h>

#include "core/Config.h"
#include "flux/FluxSource.h"

namespace nuflow {
namespace {

// Convert per-bin contents to per-GeV densities if needed.
std::vector<double> ToDensity(const std::vector<double>& edges, std::vector<double> values, bool perBin, double scale) {
  for (size_t i = 0; i < values.size(); ++i) {
    values[i] *= scale;
    if (perBin) values[i] /= (edges[i + 1] - edges[i]);
  }
  return values;
}

// ROOT TH1 (any binning). MINERvA LE 2016 file: TH1F "numu_fhc" etc., nu/m^2/1e6 POT per 0.5 GeV bin.
std::unique_ptr<HistogramFlux> ReadRoot(const std::string& file, const std::string& hist, bool perBin, double scale,
                                        const std::string& units) {
  std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "READ"));
  if (!f || f->IsZombie()) throw std::runtime_error("flux: cannot open ROOT file " + file);
  auto* h = dynamic_cast<TH1*>(f->Get(hist.c_str()));
  if (!h) throw std::runtime_error("flux: histogram '" + hist + "' not found in " + file);
  const int n = h->GetNbinsX();
  std::vector<double> edges(n + 1), values(n);
  for (int i = 1; i <= n; ++i) {
    edges[i - 1] = h->GetXaxis()->GetBinLowEdge(i);
    values[i - 1] = h->GetBinContent(i);
  }
  edges[n] = h->GetXaxis()->GetBinUpEdge(n);
  return std::make_unique<HistogramFlux>(edges, ToDensity(edges, values, perBin, scale), file + ":" + hist, units);
}

// CSV with a header line and columns BinLow,BinHigh,Flux (MINERvA 2019 nu-e constrained flux release).
std::unique_ptr<HistogramFlux> ReadCsv(const std::string& file, bool perBin, double scale, const std::string& units) {
  std::ifstream in(file);
  if (!in) throw std::runtime_error("flux: cannot open " + file);
  std::string line;
  std::vector<double> edges, values;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream ss(line);
    double lo, hi, v;
    if (!(ss >> lo >> hi >> v)) continue;  // header or malformed
    if (edges.empty()) edges.push_back(lo);
    else if (std::abs(edges.back() - lo) > 1e-9 * std::max(1.0, std::abs(lo)))
      throw std::runtime_error("flux: non-contiguous bins in " + file);
    edges.push_back(hi);
    values.push_back(v);
  }
  if (values.empty()) throw std::runtime_error("flux: no bins read from " + file);
  return std::make_unique<HistogramFlux>(edges, ToDensity(edges, values, perBin, scale), file, units);
}

// MINERvA LE 2016 text release: line 1 = flux in nu/m^2/POT per 0.5 GeV bin, 200 bins from 0 to 100 GeV;
// line 2 = uncertainties; then the covariance matrix (ignored here).
std::unique_ptr<HistogramFlux> ReadMinervaLeText(const std::string& file, double scale, const std::string& units) {
  std::ifstream in(file);
  if (!in) throw std::runtime_error("flux: cannot open " + file);
  std::string line;
  while (std::getline(in, line) && line.find_first_not_of(" \t\r") == std::string::npos) {}
  std::istringstream ss(line);
  std::vector<double> values;
  for (double v; ss >> v;) values.push_back(v);
  if (values.size() != 200) throw std::runtime_error("flux: expected 200 bins in MINERvA LE file " + file);
  std::vector<double> edges(201);
  for (int i = 0; i <= 200; ++i) edges[i] = 0.5 * i;
  return std::make_unique<HistogramFlux>(edges, ToDensity(edges, values, true, scale), file, units);
}

// Generic whitespace table: either "E_low E_high flux" or "E_center flux" (uniform bins assumed).
std::unique_ptr<HistogramFlux> ReadTable(const std::string& file, bool perBin, double scale, const std::string& units) {
  std::ifstream in(file);
  if (!in) throw std::runtime_error("flux: cannot open " + file);
  std::vector<std::vector<double>> rows;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::replace(line.begin(), line.end(), ',', ' ');
    std::istringstream ss(line);
    std::vector<double> r;
    for (double v; ss >> v;) r.push_back(v);
    if (!r.empty()) rows.push_back(r);
  }
  if (rows.size() < 2) throw std::runtime_error("flux: too few rows in " + file);
  std::vector<double> edges, values;
  if (rows[0].size() >= 3) {
    for (const auto& r : rows) {
      if (edges.empty()) edges.push_back(r[0]);
      edges.push_back(r[1]);
      values.push_back(r[2]);
    }
  } else {
    const size_t n = rows.size();
    for (size_t i = 0; i < n; ++i) {
      const double c = rows[i][0];
      const double lo = i == 0 ? c - 0.5 * (rows[1][0] - c) : 0.5 * (rows[i - 1][0] + c);
      if (i == 0) edges.push_back(lo);
      const double hi = i + 1 < n ? 0.5 * (c + rows[i + 1][0]) : c + 0.5 * (c - rows[i - 1][0]);
      edges.push_back(hi);
      values.push_back(rows[i][1]);
    }
  }
  return std::make_unique<HistogramFlux>(edges, ToDensity(edges, values, perBin, scale), file, units);
}

}  // namespace

std::unique_ptr<FluxSource> MakeFlux(const YAML::Node& node) {
  const std::string type = Get<std::string>(node, "type", "histogram");
  if (type == "mono") return std::make_unique<MonoEnergeticFlux>(node["energy"].as<double>());
  if (type == "flat") {
    const double lo = node["emin"].as<double>(), hi = node["emax"].as<double>();
    return std::make_unique<HistogramFlux>(std::vector<double>{lo, hi}, std::vector<double>{1.0},
                                           "flat " + std::to_string(lo) + "-" + std::to_string(hi) + " GeV");
  }
  if (type != "histogram") throw std::runtime_error("flux: unknown type '" + type + "'");

  const std::string format = Get<std::string>(node, "format", "root");
  const std::string file = ResolvePath(node["file"].as<std::string>());
  const std::string density = Get<std::string>(node, "density", "per_GeV");
  if (density != "per_bin" && density != "per_GeV") throw std::runtime_error("flux: density must be per_bin|per_GeV");
  const bool perBin = density == "per_bin";
  const double scale = Get<double>(node, "scale", 1.0);
  const std::string units = Get<std::string>(node, "units", "arb. units");

  std::unique_ptr<HistogramFlux> flux;
  if (format == "root") flux = ReadRoot(file, node["histogram"].as<std::string>(), perBin, scale, units);
  else if (format == "csv") flux = ReadCsv(file, perBin, scale, units);
  else if (format == "minerva_le_txt") flux = ReadMinervaLeText(file, scale, units);
  else if (format == "table") flux = ReadTable(file, perBin, scale, units);
  else throw std::runtime_error("flux: unknown format '" + format + "'");

  flux->Clip(Get<double>(node, "emin", 0.0), Get<double>(node, "emax", std::numeric_limits<double>::max()));
  return flux;
}

}  // namespace nuflow
