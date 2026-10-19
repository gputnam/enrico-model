#pragma once
// Minimal work-sharing loop over [0, n) on a fixed number of threads.

#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace nuflow {

template <typename Fn>
void ParallelFor(size_t n, int nThreads, Fn&& fn) {
  nThreads = std::max(1, std::min<int>(nThreads, static_cast<int>(n)));
  std::atomic<size_t> next{0};
  std::exception_ptr error;
  std::mutex errorMutex;
  auto worker = [&](int tid) {
    try {
      for (size_t i; (i = next.fetch_add(1)) < n;) fn(i, tid);
    } catch (...) {
      std::lock_guard<std::mutex> lock(errorMutex);
      if (!error) error = std::current_exception();
      next = n;
    }
  };
  std::vector<std::thread> threads;
  for (int t = 1; t < nThreads; ++t) threads.emplace_back(worker, t);
  worker(0);
  for (auto& t : threads) t.join();
  if (error) std::rethrow_exception(error);
}

}  // namespace nuflow
