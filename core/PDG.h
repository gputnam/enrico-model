#pragma once
// PDG Monte Carlo particle numbering helpers.

#include <cstdlib>
#include <string>

#include "core/Constants.h"

namespace nuflow::pdg {

inline constexpr int kElectron = 11, kNuE = 12, kMuon = 13, kNuMu = 14, kTau = 15, kNuTau = 16;
inline constexpr int kProton = 2212, kNeutron = 2112;

inline bool IsNeutrino(int code) {
  const int a = std::abs(code);
  return a == kNuE || a == kNuMu || a == kNuTau;
}
inline bool IsAntiNeutrino(int code) { return IsNeutrino(code) && code < 0; }

// Charged lepton produced in a CC interaction of neutrino `nu` (nu -> l-, nubar -> l+).
inline int ChargedLeptonFor(int nu) { return nu > 0 ? nu - 1 : nu + 1; }

inline double LeptonMass(int code) {
  switch (std::abs(code)) {
    case kElectron: return constants::kMassElectron;
    case kMuon: return constants::kMassMuon;
    case kTau: return constants::kMassTau;
    default: return 0.0;
  }
}

inline double NucleonMass(int code) { return code == kProton ? constants::kMassProton : constants::kMassNeutron; }

// Nuclear codes 10LZZZAAAI.
inline int IonCode(int Z, int A) {
  if (A == 1 && Z == 1) return kProton;
  if (A == 1 && Z == 0) return kNeutron;
  return 1000000000 + Z * 10000 + A * 10;
}
inline int IonZ(int code) {
  if (code == kProton) return 1;
  if (code == kNeutron) return 0;
  return (code / 10000) % 1000;
}
inline int IonA(int code) {
  if (code == kProton || code == kNeutron) return 1;
  return (code / 10) % 1000;
}

}  // namespace nuflow::pdg
