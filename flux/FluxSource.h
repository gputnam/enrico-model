#pragma once
// Neutrino flux sources. A flux is either a single energy or a binned spectrum.
// Only the spectral shape matters for event generation; the absolute normalisation is kept
// so that the flux integral can be reported (e.g. nu / m^2 / POT).

#include <memory>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

namespace nuflow {

class FluxSource {
 public:
  virtual ~FluxSource() = default;
  virtual bool IsMonoEnergetic() const = 0;
  virtual double EMin() const = 0;
  virtual double EMax() const = 0;
  virtual std::string Description() const = 0;
};

class MonoEnergeticFlux : public FluxSource {
 public:
  explicit MonoEnergeticFlux(double energy) : energy_(energy) {}
  bool IsMonoEnergetic() const override { return true; }
  double EMin() const override { return energy_; }
  double EMax() const override { return energy_; }
  double Energy() const { return energy_; }
  std::string Description() const override;

 private:
  double energy_;
};

// Piecewise-constant spectrum: bin edges [GeV] and flux density per GeV in each bin.
class HistogramFlux : public FluxSource {
 public:
  HistogramFlux(std::vector<double> edges, std::vector<double> densityPerGeV, std::string description,
                std::string units = "arb. units / GeV");
  bool IsMonoEnergetic() const override { return false; }
  double EMin() const override { return edges_.front(); }
  double EMax() const override { return edges_.back(); }
  std::string Description() const override { return description_; }
  const std::vector<double>& Edges() const { return edges_; }
  const std::vector<double>& Density() const { return density_; }
  size_t NBins() const { return density_.size(); }
  double Integral() const;  // integral over energy, in the file's units (without "/GeV")
  const std::string& Units() const { return units_; }

  // Restrict to [emin, emax] (bins partially inside are clipped).
  void Clip(double emin, double emax);

 private:
  std::vector<double> edges_;
  std::vector<double> density_;
  std::string description_;
  std::string units_;
};

// Build a flux from the `beam.flux` section of a run configuration. Supported:
//   {type: mono, energy: E}
//   {type: flat, emin: E0, emax: E1}
//   {type: histogram, format: root|csv|minerva_le_txt|table, file: ..., histogram: name,
//    density: per_bin|per_GeV, scale: s, emin: ..., emax: ...}
std::unique_ptr<FluxSource> MakeFlux(const YAML::Node& node);

}  // namespace nuflow
