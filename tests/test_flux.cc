#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "core/Config.h"
#include "flux/FluxSource.h"

using namespace nuflow;
using Catch::Approx;

namespace {
std::unique_ptr<HistogramFlux> Load(const std::string& yaml) {
  auto f = MakeFlux(YAML::Load(yaml));
  return std::unique_ptr<HistogramFlux>(dynamic_cast<HistogramFlux*>(f.release()));
}
}  // namespace

TEST_CASE("MINERvA LE 2016 flux: ROOT and text releases agree", "[flux]") {
  auto root = Load(
      "{type: histogram, format: root, file: data/flux/minerva_le_2016/minerva_flux.root, histogram: numu_fhc,"
      " density: per_bin, scale: 1.0e-6}");
  auto txt = Load("{type: histogram, format: minerva_le_txt, file: data/flux/minerva_le_2016/numuFHC.txt}");
  // Integral of the published histogram: 285.887 nu/m^2/1e6 POT (sum over 0.5 GeV bins).
  CHECK(root->Integral() == Approx(285.887e-6).epsilon(1e-4));
  CHECK(txt->Integral() == Approx(root->Integral()).epsilon(1e-5));
  REQUIRE(root->NBins() == txt->NBins());
  for (size_t i = 0; i < root->NBins(); ++i) CHECK(root->Density()[i] == Approx(txt->Density()[i]).epsilon(1e-5));
  // Per-GeV density: first bin content 2.57268 per 0.5 GeV per 1e6 POT.
  CHECK(root->Density()[0] == Approx(2.57268e-6 / 0.5).epsilon(1e-5));
}

TEST_CASE("MINERvA ME 2019 CSV flux", "[flux]") {
  auto csv = Load("{type: histogram, format: csv, file: data/flux/minerva_me_2019/MINERvA_Flux_pdg14_500MeVBins.csv}");
  CHECK(csv->EMin() == Approx(0.0));
  CHECK(csv->Edges()[1] == Approx(0.5));
  CHECK(csv->Density()[0] == Approx(4.501604e-06));
  CHECK(csv->EMax() <= 100.0);
}

TEST_CASE("Flux clipping and mono-energetic beams", "[flux]") {
  auto f = Load(
      "{type: histogram, format: root, file: data/flux/minerva_le_2016/minerva_flux.root, histogram: numu_fhc,"
      " density: per_bin, emin: 1.25, emax: 10.0}");
  CHECK(f->EMin() == Approx(1.25));
  CHECK(f->EMax() == Approx(10.0));
  CHECK(f->Edges()[1] == Approx(1.5));
  auto mono = MakeFlux(YAML::Load("{type: mono, energy: 2.5}"));
  CHECK(mono->IsMonoEnergetic());
  CHECK(mono->EMin() == 2.5);
}
