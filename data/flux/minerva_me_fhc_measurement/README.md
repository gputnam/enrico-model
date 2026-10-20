# MINERvA ME FHC measurement flux (numu)

`flux_pdg14.csv` is the flux that the measurements in `measurement/` divide by: the MINERvA
ME FHC CV flux (nu-e scattering constrained, `MATFluxAndReweightFiles`), weighted by the data
POT of the 12 ME FHC playlists 1A-1P of the v1 dataframes, 0-100 GeV (128 bins).
Integral: 6.24334e-04 nu/m^2/POT, identical to `flux_integral` of every measurement result.

Produced with (in measurement/):

    venv/bin/python -m xsec.export_flux -o flux_pdg14.csv --tune v1

Use this file (not `minerva_me_2019/`) whenever a prediction is compared with measurement/.
