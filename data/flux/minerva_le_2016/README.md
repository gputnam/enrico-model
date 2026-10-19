# MINERvA low-energy (LE) NuMI flux, 2016

**Source:** L. Aliaga *et al.* (MINERvA), "Neutrino Flux Predictions for the NuMI Beam", *Phys. Rev. D* **94**, 092005 (2016), [arXiv:1607.00704](https://arxiv.org/abs/1607.00704).
**Files:** downloaded unmodified from the arXiv ancillary files, https://arxiv.org/src/1607.00704v2/anc/.

| File | Contents |
|---|---|
| `minerva_flux.root` | `TH1F` `{numu,numubar,nue,nuebar}_{fhc,rhc}`: 200 bins, 0–100 GeV, in **ν / m² / 10⁶ POT per 0.5 GeV bin**; also `*_covmx` and `*_cormx` `TH2F`s |
| `{numu,numubar,nue,nuebar}{FHC,RHC}.txt` | Line 1: flux in ν / m² / POT per 0.5 GeV bin (200 bins). Line 2: uncertainty. Then the 200×200 covariance |
| `read_flux.py` | The collaboration's reader script, which produced the ROOT file from the text files |

**Use in nuflow:** `format: root, histogram: numu_fhc, density: per_bin, scale: 1.0e-6`. That gives ν / m² / POT / GeV. The integrated `numu_fhc` flux is 2.859×10⁻⁴ ν / m² / POT.
