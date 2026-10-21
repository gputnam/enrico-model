#include "kinematics/MECKinematics.h"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "core/Constants.h"
#include "core/PDG.h"

namespace nuflow {

using namespace constants;
using cplx = std::complex<double>;

std::string MECKinematics::Name() const {
  std::ostringstream s;
  s.precision(6);
  s << "CCMEC[" << current_.Name() << ",qscale=" << opt_.proposalQScale << ",q3max=" << opt_.q3Max << "]";
  return s.str();
}

double MECKinematics::Generate(int nuPdg, double enu, const Nucleus& nucleus, const InitialState& initial,
                               const NuclearModel& nuclearModel, Event& event, Random& rng) const {
  if (initial.nHit != 2) return 0;
  const bool nu = nuPdg > 0;
  const int helicity = nu ? +1 : -1;
  const int leptonPdg = pdg::ChargedLeptonFor(nuPdg);
  const double mLep = pdg::LeptonMass(leptonPdg);
  const double M = kMassNucleon;

  // Charge channel, hole labelling (h1, h2) and final species (f1, f2); see MECCurrent.h.
  const HitNucleon* hA = &initial.hit[0];
  const HitNucleon* hB = &initial.hit[1];
  const int nProtons = (hA->pdg == pdg::kProton) + (hB->pdg == pdg::kProton);
  MECChannel ch;
  int pdg1, pdg2;
  double kFpair;  // Fermi momentum entering the pair density (theory/mec.md, Eq. for P_X)
  if (nu) {
    if (nProtons == 1) {
      ch = MECChannel::kNuNPtoPP;
      if (hA->pdg != pdg::kNeutron) std::swap(hA, hB);  // h1 = neutron
      pdg1 = pdg2 = pdg::kProton;
      kFpair = nucleus.KFermi(pdg::kProton);
    } else if (nProtons == 0) {
      ch = MECChannel::kNuNNtoNP;
      pdg1 = pdg::kNeutron;
      pdg2 = pdg::kProton;
      kFpair = nucleus.KFermi(pdg::kNeutron);
    } else {
      return 0;
    }
  } else {
    if (nProtons == 2) {
      ch = MECChannel::kNubarPPtoNP;
      pdg1 = pdg::kNeutron;
      pdg2 = pdg::kProton;
      kFpair = nucleus.KFermi(pdg::kProton);
    } else if (nProtons == 1) {
      ch = MECChannel::kNubarNPtoNN;
      if (hA->pdg != pdg::kNeutron) std::swap(hA, hB);
      pdg1 = pdg2 = pdg::kNeutron;
      kFpair = nucleus.KFermi(pdg::kNeutron);
    } else {
      return 0;
    }
  }
  const double m1 = pdg::NucleonMass(pdg1), m2 = pdg::NucleonMass(pdg2);

  // ---- Lepton: proposal in (|q|, omega, phi); d3k'/E' = (|q|/E_nu) d omega d|q| d phi.
  const double eps = enu;
  if (eps <= mLep * 1.0001) return 0;
  const double kmax = std::sqrt(eps * eps - mLep * mLep);
  const double qlo = eps - kmax;
  double qhi = eps + kmax;
  if (opt_.q3Max > 0) qhi = std::min(qhi, opt_.q3Max);
  if (qhi <= qlo) return 0;
  const double q0s = opt_.proposalQScale;
  const double a = 1 / (qlo + q0s), b = 1 / (qhi + q0s);
  const double q = 1 / (a - rng.Uniform() * (a - b)) - q0s;
  const double gq = 1 / ((q + q0s) * (q + q0s) * (a - b));

  const FourVector pairIn = hA->Momentum() + hB->Momentum();  // off-shell (removal energies included)
  const double Ppair = pairIn.P();
  const double wHi = eps - std::sqrt((eps - q) * (eps - q) + mLep * mLep);
  const double qMinPair = std::max(0.0, q - Ppair);
  const double wLo = std::max(0.0, std::sqrt((m1 + m2) * (m1 + m2) + qMinPair * qMinPair) - pairIn.e);
  if (wHi <= wLo) return 0;
  const double w = wLo + (wHi - wLo) * rng.Uniform();
  const double g = gq / (wHi - wLo);

  const double eLep = eps - w;
  const double kp = std::sqrt(std::max(0.0, eLep * eLep - mLep * mLep));
  if (kp <= 0) return 0;
  const double cosT = std::clamp((eps * eps + kp * kp - q * q) / (2 * eps * kp), -1.0, 1.0);
  const double sinT = std::sqrt(std::max(0.0, 1 - cosT * cosT));
  const double phi = 2 * kPi * rng.Uniform();
  const FourVector k(eps, 0, 0, eps);
  const FourVector lepton(eLep, kp * sinT * std::cos(phi), kp * sinT * std::sin(phi), kp * cosT);
  const FourVector qv = k - lepton;

  // ---- Outgoing pair: two-body decay of P = q + h1 + h2 in its rest frame (isotropic proposal).
  const FourVector P = qv + pairIn;
  const double s = P.M2();
  if (P.e <= 0 || s <= (m1 + m2) * (m1 + m2)) return 0;
  const double W = std::sqrt(s);
  const double lam = (s - (m1 + m2) * (m1 + m2)) * (s - (m1 - m2) * (m1 - m2));
  const double pstar = std::sqrt(std::max(0.0, lam)) / (2 * W);
  const double c = 2 * rng.Uniform() - 1, sn = std::sqrt(std::max(0.0, 1 - c * c)), ph = 2 * kPi * rng.Uniform();
  const FourVector p1cm(std::sqrt(pstar * pstar + m1 * m1), pstar * sn * std::cos(ph), pstar * sn * std::sin(ph),
                        pstar * c);
  const FourVector p1 = p1cm.Boosted(P.BoostVector());
  const FourVector p2 = P - p1;
  if (nuclearModel.PauliBlocked(nucleus, pdg1, p1.Vect()) || nuclearModel.PauliBlocked(nucleus, pdg2, p2.Vect()))
    return 0;

  // ---- Matrix element: on-shell nucleons of mass M, effective transfer Q~ = p1' + p2' - h1 - h2
  //      (energy transfer reduced by the pair removal energy, the de Forest-type prescription of QE).
  auto onShell = [&](const ThreeVector& p) { return FourVector(std::sqrt(p.Mag2() + M * M), p); };
  const FourVector h1 = onShell(hA->p), h2 = onShell(hB->p);
  const FourVector f1 = onShell(p1.Vect()), f2 = onShell(p2.Vect());
  const FourVector Qt = f1 + f2 - h1 - h2;
  MECCurrent::Amplitudes A;
  current_.Antisymmetrized(ch, Qt, h1, h2, f1, f2, A);
  cplx L[4][4];
  LeptonTensor(k, lepton, mLep, helicity, L);
  double S = 0;
  for (int i = 0; i < 16; ++i) S += ContractLepton(L, A.a[i]);
  if (S <= 0) return 0;

  // sigma = P_X < G^2 Vud^2 q p* S / (256 pi^3 E_nu^2 E1 E2 W g) >, P_X = kF^3/(12 pi^2) per unit multiplicity.
  const double pairDensity = kFpair * kFpair * kFpair / (12 * kPi * kPi);
  const double weight = pairDensity * kGF * kGF * kVud * kVud * q * pstar * S /
                        (256 * kPi * kPi * kPi * eps * eps * h1.e * h2.e * W * g) * kGeV2To1e38Cm2;

  // ---- Event record.
  event.nuPdg = nuPdg;
  event.targetPdg = nucleus.pdg;
  event.enu = enu;
  event.Q2 = -qv.M2();
  event.q0 = qv.e;
  event.q3 = qv.P();
  event.W = std::sqrt(std::max(0.0, M * M + 2 * M * qv.e - event.Q2));  // free-nucleon-at-rest convention
  event.initial = initial;
  event.particles.clear();
  event.Add(nuPdg, kStatusBeam, k);
  event.Add(nucleus.pdg, kStatusTarget, FourVector(nucleus.mass, 0, 0, 0));
  const int i1 = event.Add(hA->pdg, kStatusStruckNucleon, hA->Momentum(), 1);
  const int i2 = event.Add(hB->pdg, kStatusStruckNucleon, hB->Momentum(), 1);
  event.Add(leptonPdg, kStatusFinal, lepton, 0);
  event.Add(pdg1, kStatusPreFSI, p1, i1);
  event.Add(pdg2, kStatusPreFSI, p2, i2);
  if (initial.remnantPdg != 0) event.Add(initial.remnantPdg, kStatusRemnant, initial.remnant, 1);
  return weight;
}

}  // namespace nuflow
