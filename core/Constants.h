#pragma once
// Physical constants and particle masses. Natural units: GeV, c = hbar = 1.

namespace nuflow::constants {

inline constexpr double kPi = 3.14159265358979323846;

// Weak interaction (PDG 2024)
inline constexpr double kGF = 1.1663788e-5;       // Fermi constant [GeV^-2]
inline constexpr double kVud = 0.97367;           // |V_ud|

// Conversion: (hbar c)^2 = 0.3893794 GeV^2 mb  ->  1 GeV^-2 = 0.3893794e-27 cm^2
inline constexpr double kGeV2ToCm2 = 0.38937937e-27;
inline constexpr double kGeV2To1e38Cm2 = kGeV2ToCm2 / 1e-38;   // GeV^-2 -> 1e-38 cm^2

// Masses [GeV] (PDG 2024)
inline constexpr double kMassProton = 0.93827209;
inline constexpr double kMassNeutron = 0.93956542;
inline constexpr double kMassNucleon = 0.5 * (kMassProton + kMassNeutron);
inline constexpr double kMassElectron = 0.51099895e-3;
inline constexpr double kMassMuon = 0.1056583755;
inline constexpr double kMassTau = 1.77686;
inline constexpr double kMassPionCharged = 0.13957039;
inline constexpr double kAtomicMassUnit = 0.93149410;

// Nucleon magnetic moments [nuclear magnetons]
inline constexpr double kMuProton = 2.79284734;
inline constexpr double kMuNeutron = -1.91304276;

}  // namespace nuflow::constants
