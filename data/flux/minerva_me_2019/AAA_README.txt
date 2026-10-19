This file describes how to use the auxiliary files in this paper to access the data contained in it.

There are two types of contained data: (1) the neutrino-electron scattering result, and (2) the neutrino flux after the constraint.

1: Neutrino electron scattering result

This is in the file
nuEscatDataReleaseFile.csv

It includes the results for the number of neutrino-electron scattering events in the six electron energy bins, the number of protons on target that generated the flux which resulted in the neutrino-electron scattering events, number of electron targets in the fiducial volume of the measurement, and fiducial area over which the flux is normalized.

2. NuMI FHC Flux

The remainder of the files are the resulting neutrino-electron scattering constrained FHC flux observed in the fiducial volume of the MINERvA scintillator tracker.

Units for the flux are neutrinos/m^2/proton on target(POT)/GeV.

Note the number of POT may be different in different analyses.  The one for this analysis of neutrino-electron scattering is given in the nuEscatDataReleaseFile.csv file.

The flux is reported in four different neutrino species, referred to by their PDG codes in the file name:
"pdg-12" electron anti-neutrinos
"pdg12" electron neutrinos
"pdg-14" muon anti-neutrinos
"pdg14" muon neutrinos

The primary files are the .csv files, which contain comma separated value representations of the data.  Start there.

There are two different types of files:
"MINERvA_Flux_pdg*" files contain the flux themselves, with column headings defining each 500 MeV wide bin of flux and the content.
"MINERvA_Flux_covariance_pdg*" contain the covariance matrices that describe the flux uncertainties.

The .root files contain root histograms with the same information in the similarly named .csv files.  The .png files are graphical representations.  