#pragma once
// Minimal, allocation-free Minkowski four-vector (metric +,-,-,-).
// Deliberately lighter than TLorentzVector: everything is inline and constexpr-friendly.

#include <cmath>

namespace nuflow {

struct ThreeVector {
  double x = 0, y = 0, z = 0;

  constexpr ThreeVector() = default;
  constexpr ThreeVector(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  constexpr ThreeVector operator+(const ThreeVector& o) const { return {x + o.x, y + o.y, z + o.z}; }
  constexpr ThreeVector operator-(const ThreeVector& o) const { return {x - o.x, y - o.y, z - o.z}; }
  constexpr ThreeVector operator*(double s) const { return {x * s, y * s, z * s}; }
  constexpr ThreeVector operator-() const { return {-x, -y, -z}; }
  constexpr double Dot(const ThreeVector& o) const { return x * o.x + y * o.y + z * o.z; }
  constexpr ThreeVector Cross(const ThreeVector& o) const {
    return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
  }
  constexpr double Mag2() const { return x * x + y * y + z * z; }
  double Mag() const { return std::sqrt(Mag2()); }
  ThreeVector Unit() const {
    const double m = Mag();
    return m > 0 ? (*this) * (1.0 / m) : ThreeVector{0, 0, 1};
  }
};

struct FourVector {
  double e = 0, px = 0, py = 0, pz = 0;

  constexpr FourVector() = default;
  constexpr FourVector(double e_, double px_, double py_, double pz_) : e(e_), px(px_), py(py_), pz(pz_) {}
  constexpr FourVector(double e_, const ThreeVector& p) : e(e_), px(p.x), py(p.y), pz(p.z) {}

  constexpr ThreeVector Vect() const { return {px, py, pz}; }
  constexpr FourVector operator+(const FourVector& o) const { return {e + o.e, px + o.px, py + o.py, pz + o.pz}; }
  constexpr FourVector operator-(const FourVector& o) const { return {e - o.e, px - o.px, py - o.py, pz - o.pz}; }
  constexpr FourVector operator*(double s) const { return {e * s, px * s, py * s, pz * s}; }
  constexpr double Dot(const FourVector& o) const { return e * o.e - px * o.px - py * o.py - pz * o.pz; }
  constexpr double M2() const { return Dot(*this); }
  double M() const {
    const double m2 = M2();
    return m2 >= 0 ? std::sqrt(m2) : -std::sqrt(-m2);
  }
  constexpr double P2() const { return px * px + py * py + pz * pz; }
  double P() const { return std::sqrt(P2()); }
  // Velocity of the rest frame of this four-vector (requires a time-like vector).
  constexpr ThreeVector BoostVector() const { return {px / e, py / e, pz / e}; }

  // Active boost by velocity beta: returns this vector as seen in a frame moving with -beta.
  FourVector Boosted(const ThreeVector& b) const {
    const double b2 = b.Mag2();
    if (b2 <= 0) return *this;
    const double gamma = 1.0 / std::sqrt(1.0 - b2);
    const double bp = b.x * px + b.y * py + b.z * pz;
    const double g2 = (gamma - 1.0) / b2;
    const double f = g2 * bp + gamma * e;
    return {gamma * (e + bp), px + f * b.x, py + f * b.y, pz + f * b.z};
  }
};

inline constexpr double Dot(const FourVector& a, const FourVector& b) { return a.Dot(b); }

}  // namespace nuflow
