#include "kinematics/BodekYang.h"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "core/Constants.h"
#include "kinematics/GRV98LO.h"

namespace nuflow {

using namespace constants;

std::string BodekYang::Name() const {
  std::ostringstream s;
  s.precision(6);
  s << "BY2013(A=" << par_.A << ",B=" << par_.B << ",N=" << par_.N << (par_.axialTypeII ? ",II" : ",I")
    << (par_.charm ? ",c" : "") << (par_.duCorrection ? ",du" : "") << ")";
  return s.str();
}

double BodekYang::XiW(double x, double Q2, double M, double mf) const {
  const double num = 2 * x * (Q2 + mf * mf + par_.B);
  const double den = Q2 * (1 + std::sqrt(1 + 4 * M * M * x * x / Q2)) + 2 * par_.A * x;
  return num / den;
}

double BodekYang::R1998(double x, double Q2) {
  // E143, Phys. Lett. B 452, 194 (1999), Eqs. 2-3 and Table II; R1998 = (Ra + Rb + Rc)/3.
  static constexpr double a[6] = {0.0485, 0.5470, 2.0621, -0.3804, 0.5090, -0.0285};
  static constexpr double b[6] = {0.0481, 0.6114, -0.3509, -0.4611, 0.7172, -0.0317};
  static constexpr double c[6] = {0.0577, 0.4644, 1.8288, 12.3708, -43.1043, 41.7415};
  x = std::max(x, 1e-4);
  const double Q4 = Q2 * Q2;
  const double theta = 1 + 12 * Q2 / (Q2 + 1) * 0.125 * 0.125 / (0.125 * 0.125 + x * x);
  const double lnQ = std::log(Q2 / 0.04);
  const double Ra = a[0] / lnQ * theta + a[1] / std::pow(Q4 * Q4 + std::pow(a[2], 4), 0.25) *
                                             (1 + a[3] * x + a[4] * x * x) * std::pow(x, a[5]);
  const double Rb = b[0] / lnQ * theta + (b[1] / Q2 + b[2] / (Q4 + 0.3 * 0.3)) * (1 + b[3] * x + b[4] * x * x) *
                                             std::pow(x, b[5]);
  const double q2thr = c[3] * x + c[4] * x * x + c[5] * x * x * x;
  const double Rc = c[0] / lnQ * theta + c[1] / std::sqrt((Q2 - q2thr) * (Q2 - q2thr) + c[2] * c[2]);
  return (Ra + Rb + Rc) / 3;
}

double BodekYang::R(double x, double Q2) {
  // Bodek-Yang Sec. 7: R1998 for Q2 > 0.3 GeV^2; below, 3.633 Q2/(Q4 + 1) R1998(x, 0.3) -> 0 at Q2 = 0.
  if (Q2 >= 0.3) return std::max(0.0, R1998(x, Q2));
  return std::max(0.0, 3.633 * Q2 / (Q2 * Q2 + 1) * R1998(x, 0.3));
}

namespace {

// Effective-LO parton content of the struck nucleon at one xi, with the K factors attached.
struct Flavours {
  double uv, dv, us, ds, s;  // xi * q(xi): valence u, d; sea u (= ubar), sea d (= dbar), s (= sbar)
};

Flavours Partons(double xi, double Q2, bool proton, bool du, double N) {
  Flavours f{0, 0, 0, 0, 0};
  if (!(xi > 0) || xi >= 1) return f;
  const PartonDensities p = GRV98LO::Instance().Evaluate(xi, Q2);
  double uv = p.uv, dv = p.dv;
  if (du && uv + dv > 0) {
    // NMC d/u correction (Bodek-Yang Eqs. 35-37): (d/u)' = d/u + delta, u_v + d_v unchanged.
    const double delta = -0.00817 + 0.0506 * xi + 0.0798 * xi * xi;
    const double den = 1 + delta * uv / (uv + dv);
    const double uvp = uv / den, dvp = (dv + uv * delta) / den;
    uv = uvp;
    dv = dvp;
  }
  f = {uv, dv, p.ubar, p.dbar, p.s};
  if (!proton) {  // isospin: u <-> d
    std::swap(f.uv, f.dv);
    std::swap(f.us, f.ds);
  }
  f.uv *= N;
  f.dv *= N;
  f.us *= N;
  f.ds *= N;
  f.s *= N;
  return f;
}

}  // namespace

StructureFunctions BodekYang::Compute(bool neutrino, bool proton, double x, double Q2, double nu, double M) const {
  StructureFunctions out;
  if (!(x > 0) || x >= 1 || !(Q2 > 0) || !(nu > 0)) return out;
  const auto& P = par_;

  // K factors (Q2 dependence only). They belong to the parton distribution: for the neutron (isospin
  // rotated proton PDFs) the "u" factors go with the distribution that was u in the proton.
  const double GD = 1 / ((1 + Q2 / P.MDipoleV2) * (1 + Q2 / P.MDipoleV2));
  const double KLW = (nu * nu + P.ClowNu) / (nu * nu);
  double KvU = KLW * (1 - GD * GD) * (Q2 + P.Cv2u) / (Q2 + P.Cv1u);
  double KvD = KLW * (1 - GD * GD) * (Q2 + P.Cv2d) / (Q2 + P.Cv1d);
  double KsU = Q2 / (Q2 + P.CseaU), KsD = Q2 / (Q2 + P.CseaD);
  const double KsS = Q2 / (Q2 + P.CseaS);
  if (!proton) {
    std::swap(KvU, KvD);
    std::swap(KsU, KsD);
  }
  double KvA = (Q2 + P.PvalAxial) / (Q2 + P.CvalAxial), KsA = (Q2 + P.PseaAxial) / (Q2 + P.CseaAxial);

  const double Q2pdf = std::max(Q2, GRV98LO::kQ2Min);
  const double tm = 1 + 4 * M * M * x * x / Q2;
  const double Vud2 = kVud * kVud, Vus2 = kVus * kVus, Vcd2 = kVcd * kVcd, Vcs2 = kVcs * kVcs;

  // Accumulates  sum K^V xi q,  sum K^A xi q  for quarks and antiquarks separately.
  struct Sum {
    double qV = 0, qA = 0, qbV = 0, qbA = 0;
  };
  auto add = [&](Sum& s, bool anti, double xq, double w, double KV, double KAx) {
    const double KA = P.axialTypeII ? KAx : KV;
    if (anti) {
      s.qbV += w * KV * xq;
      s.qbA += w * KA * xq;
    } else {
      s.qV += w * KV * xq;
      s.qA += w * KA * xq;
    }
  };

  // --- light final quark (M_f = 0)
  const double xi0 = XiW(x, Q2, M, 0.0);
  Sum light;
  if (xi0 < 1) {
    const Flavours f = Partons(xi0, Q2pdf, proton, P.duCorrection, P.N);
    if (neutrino) {  // W+ : d, s -> u ; ubar -> dbar, sbar
      add(light, false, f.dv, Vud2, KvD, KvA);
      add(light, false, f.ds, Vud2, KsD, KsA);
      add(light, false, f.s, Vus2, KsS, KsA);
      add(light, true, f.us, Vud2 + Vus2, KsU, KsA);
    } else {  // W- : u -> d, s ; dbar, sbar -> ubar
      add(light, false, f.uv, Vud2 + Vus2, KvU, KvA);
      add(light, false, f.us, Vud2 + Vus2, KsU, KsA);
      add(light, true, f.ds, Vud2, KsD, KsA);
      add(light, true, f.s, Vus2, KsS, KsA);
    }
  }
  // --- charm final quark (slow rescaling, M_f = M_c)
  Sum heavy;
  double Kc = 0, xic = 1;
  if (P.charm) {
    xic = XiW(x, Q2, M, P.charmMass);
    Kc = Q2 / (Q2 + P.charmMass * P.charmMass);
    if (xic < 1) {
      const Flavours f = Partons(xic, Q2pdf, proton, P.duCorrection, P.N);
      if (neutrino) {  // d, s -> c
        add(heavy, false, f.dv, Vcd2, KvD, KvA);
        add(heavy, false, f.ds, Vcd2, KsD, KsA);
        add(heavy, false, f.s, Vcs2, KsS, KsA);
      } else {  // dbar, sbar -> cbar
        add(heavy, true, f.ds, Vcd2, KsD, KsA);
        add(heavy, true, f.s, Vcs2, KsS, KsA);
      }
    }
  }

  const double H = BodekYang::H(x);
  // F2 = F2^V + F2^A; 2xF1 = 2 (2xF1^V) with 2xF1^A = 2xF1^V (Sec. 11); xF3 = 2 H (sum q - sum qbar).
  const double F2Vl = light.qV + light.qbV, F2Al = light.qA + light.qbA;
  const double F2Vh = heavy.qV + heavy.qbV, F2Ah = heavy.qA + heavy.qbA;
  out.F2 = F2Vl + F2Al + F2Vh + F2Ah;
  double twoxF1 = 2 * F2Vl * tm / (1 + R(x, Q2));
  if (F2Vh > 0) twoxF1 += 2 * Kc * F2Vh * tm / (1 + R(xic, Q2));
  const double xF3 = 2 * H * (light.qV - light.qbV) + 2 * H * Kc * (heavy.qV - heavy.qbV);
  out.F1 = twoxF1 / (2 * x);
  out.F3 = xF3 / x;
  out.F4 = 0;
  out.F5 = out.F2 / (2 * x);
  return out;
}

}  // namespace nuflow
