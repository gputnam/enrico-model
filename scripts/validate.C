// Validation plots for a nuflow output file.
//   root -l -b -q 'scripts/validate.C("runs/minerva_le_fhc_numu_CH.root", "runs/plots")'
// Produces PNGs and prints summary statistics, including a check of four-momentum conservation and
// of the event energy spectrum against flux x sigma.

#include <cmath>
#include <cstdio>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

void validate(const char* path = "runs/minerva_le_fhc_numu_CH.root", const char* outdir = "runs/plots") {
  gStyle->SetOptStat(0);
  gSystem->mkdir(outdir, true);
  TFile f(path);
  auto* t = static_cast<TTree*>(f.Get("events"));
  auto* run = static_cast<TTree*>(f.Get("run"));
  double sigmaAvg = 0;
  run->SetBranchAddress("flux_avg_xsec", &sigmaAvg);
  run->GetEntry(0);

  const int kMax = 64;
  int n, pdg[kMax], status[kMax], target, nu;
  double E[kMax], px[kMax], py[kMax], pz[kMax], enu, Q2, hitp[3];
  t->SetBranchAddress("n", &n);
  t->SetBranchAddress("pdg", pdg);
  t->SetBranchAddress("status", status);
  t->SetBranchAddress("E", E);
  t->SetBranchAddress("px", px);
  t->SetBranchAddress("py", py);
  t->SetBranchAddress("pz", pz);
  t->SetBranchAddress("Enu", &enu);
  t->SetBranchAddress("Q2", &Q2);
  t->SetBranchAddress("target_pdg", &target);
  t->SetBranchAddress("nu_pdg", &nu);
  t->SetBranchAddress("hit_p", hitp);

  TH1D hE("hE", ";E_{#nu} [GeV];events", 80, 0, 20);
  TH1D hQ2("hQ2", ";Q^{2} [GeV^{2}];events", 100, 0, 2.5);
  TH1D hPmu("hPmu", ";p_{#mu} [GeV];events", 100, 0, 20);
  TH1D hCos("hCos", ";cos#theta_{#mu};events", 100, -1, 1);
  TH1D hTp("hTp", ";nucleon kinetic energy [GeV];events", 100, 0, 2.5);
  TH1D hKF("hKF", ";struck-nucleon |p| [GeV];events", 100, 0, 0.3);
  long long nEv = t->GetEntries(), nH = 0;
  double maxViolation = 0;
  for (long long i = 0; i < nEv; ++i) {
    t->GetEntry(i);
    hE.Fill(enu);
    hQ2.Fill(Q2);
    if (target == 2212) ++nH;
    else hKF.Fill(std::sqrt(hitp[0] * hitp[0] + hitp[1] * hitp[1] + hitp[2] * hitp[2]));
    double in[4] = {0, 0, 0, 0}, out[4] = {0, 0, 0, 0};
    for (int j = 0; j < n; ++j) {
      double* acc = (status[j] == 4 || status[j] == 20) ? in : (status[j] == 1 || status[j] == 23) ? out : nullptr;
      if (acc) acc[0] += E[j], acc[1] += px[j], acc[2] += py[j], acc[3] += pz[j];
      if (status[j] == 1 && std::abs(pdg[j]) == 13) {
        const double p = std::sqrt(px[j] * px[j] + py[j] * py[j] + pz[j] * pz[j]);
        hPmu.Fill(p);
        hCos.Fill(pz[j] / p);
      }
      if (status[j] == 1 && (pdg[j] == 2212 || pdg[j] == 2112))
        hTp.Fill(E[j] - std::sqrt(E[j] * E[j] - px[j] * px[j] - py[j] * py[j] - pz[j] * pz[j]));
    }
    for (int k = 0; k < 4; ++k) maxViolation = std::max(maxViolation, std::abs(in[k] - out[k]));
  }

  // Expected E_nu shape: flux x sigma_tot (sum over xsec graphs weighted by abundance is already
  // folded into the generator; here we use the graphs and the target abundances from the file).
  auto* flux = static_cast<TH1D*>(f.Get("flux"));
  TH1D hExp("hExp", "", 80, 0, 20);
  if (flux) {
    TList* keys = f.GetListOfKeys();
    for (int b = 1; b <= hExp.GetNbinsX(); ++b) {
      double s = 0;
      const double lo = hExp.GetXaxis()->GetBinLowEdge(b), hi = hExp.GetXaxis()->GetBinUpEdge(b);
      for (int j = 0; j < 50; ++j) {
        const double e = lo + (hi - lo) * (j + 0.5) / 50;
        double sig = 0;
        for (TObject* k : *keys) {
          TString name = k->GetName();
          if (!name.BeginsWith("xsec_")) continue;
          sig += static_cast<TGraph*>(f.Get(name))->Eval(e);  // counts are 1 per nucleus in CH
        }
        const int fb = flux->FindBin(e);
        if (fb >= 1 && fb <= flux->GetNbinsX()) s += flux->GetBinContent(fb) * sig;
      }
      hExp.SetBinContent(b, s);
    }
    hExp.Scale(hE.Integral() / hExp.Integral());
  }

  std::printf("events: %lld, on hydrogen: %lld (%.2f%%)\n", nEv, nH, 100.0 * nH / nEv);
  std::printf("flux-averaged sigma: %.4f x 1e-38 cm^2 per target unit\n", sigmaAvg);
  std::printf("max |sum p_in - sum p_out| over events: %.3g GeV\n", maxViolation);
  std::printf("<E_nu> = %.3f GeV, <Q2> = %.3f GeV^2, <p_mu> = %.3f GeV\n", hE.GetMean(), hQ2.GetMean(), hPmu.GetMean());
  if (flux) {
    double chi2 = 0;
    int ndf = 0;
    for (int b = 1; b <= hE.GetNbinsX(); ++b) {
      const double e = hExp.GetBinContent(b);
      if (e < 20) continue;
      chi2 += std::pow(hE.GetBinContent(b) - e, 2) / e;
      ++ndf;
    }
    std::printf("E_nu shape vs flux x sigma: chi2/ndf = %.1f/%d\n", chi2, ndf);
  }

  TCanvas c("c", "", 1500, 900);
  c.Divide(3, 2);
  c.cd(1);
  hE.SetLineWidth(2);
  hE.Draw("hist");
  if (flux) {
    hExp.SetLineColor(kRed);
    hExp.SetLineStyle(2);
    hExp.Draw("hist same");
    TLegend* l = new TLegend(0.5, 0.75, 0.88, 0.88);
    l->AddEntry(&hE, "generated", "l");
    l->AddEntry(&hExp, "#Phi(E) #sigma(E)", "l");
    l->Draw();
  }
  TH1D* hs[] = {&hQ2, &hPmu, &hCos, &hTp, &hKF};
  for (int i = 0; i < 5; ++i) {
    c.cd(i + 2);
    hs[i]->SetLineWidth(2);
    hs[i]->Draw("hist");
  }
  c.SaveAs(Form("%s/validation.png", outdir));

  TCanvas c2("c2", "", 800, 600);
  bool first = true;
  int color = 1;
  for (TObject* k : *f.GetListOfKeys()) {
    TString name = k->GetName();
    if (!name.BeginsWith("xsec_")) continue;
    auto* g = static_cast<TGraph*>(f.Get(name));
    g->SetLineColor(color++);
    g->SetLineWidth(2);
    g->Draw(first ? "AL" : "L");
    first = false;
  }
  c2.SetLogx();
  c2.SaveAs(Form("%s/xsec.png", outdir));
}
