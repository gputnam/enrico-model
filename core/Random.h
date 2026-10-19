#pragma once
// xoshiro256++ pseudo-random generator (Blackman & Vigna, 2019), seeded through SplitMix64.
// Each event batch / cross-section-table node gets its own stream derived from (seed, stream id),
// which makes results independent of the number of worker threads.

#include <cmath>
#include <cstdint>

namespace nuflow {

class Random {
 public:
  explicit Random(uint64_t seed = 0x9E3779B97F4A7C15ULL, uint64_t stream = 0) { Seed(seed, stream); }

  void Seed(uint64_t seed, uint64_t stream = 0) {
    uint64_t sm = seed ^ (0xD1B54A32D192ED03ULL * (stream + 1));
    for (auto& w : s_) w = SplitMix64(sm);
  }

  uint64_t Next() {
    const uint64_t result = Rotl(s_[0] + s_[3], 23) + s_[0];
    const uint64_t t = s_[1] << 17;
    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = Rotl(s_[3], 45);
    return result;
  }

  // Uniform in [0, 1).
  double Uniform() { return (Next() >> 11) * 0x1.0p-53; }
  // Uniform in (0, 1]; safe for log().
  double UniformOpen() { return ((Next() >> 11) + 1) * 0x1.0p-53; }
  double Uniform(double lo, double hi) { return lo + (hi - lo) * Uniform(); }

 private:
  static uint64_t Rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  static uint64_t SplitMix64(uint64_t& x) {
    uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
  }
  uint64_t s_[4];
};

}  // namespace nuflow
