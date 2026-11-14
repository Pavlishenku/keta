// Persistent worker pool for fine-grained, deterministic loops.
//
//   parallel_for(n, threads, [&](std::size_t i) { out[i] = f(i); });
//
// runs body(i) for every i in [0, n), each index exactly once, on up to
// `threads` workers (0 = hardware concurrency, 1 = inline). Determinism is the
// caller's contract and is kept by construction when every index writes only
// its own pre-allocated slot and any reduction is done afterwards, in index
// order, on the calling thread: the thread count and the schedule then cannot
// change a single ULP.
//
// The workers are started once and reused (a staged-construction run makes
// tens of thousands of small parallel calls; spawning threads per call would
// cost more than the work). A call made from inside a worker runs inline, so
// nested use cannot deadlock. The first exception (lowest index) is rethrown
// on the calling thread after every index has run.
#pragma once
#include <cstddef>
#include <functional>

namespace keta {

void parallel_for(std::size_t n, int threads, const std::function<void(std::size_t)>& body);

// The worker count `threads` resolves to (0 = hardware concurrency).
int resolve_threads(int threads);

}  // namespace keta
