// Cooperative cancellation and progress hooks for the heavy solve drivers.
//
// A driver whose options carry a non-null `hooks` polls the cancel token at
// phase / iteration boundaries only — never inside a factorization or an
// eigensolver kernel — and throws `Cancelled` when it fired; between polls it
// reports a progress fraction in [0, 1] with a short phase label. Without
// hooks (the default, a null pointer) every driver behaves exactly as before:
// the numbers are bit-identical, polling only reads a flag.
//
// Polling points, per driver:
//   static (cached)      entry, "assembly" 0, "factorization" 0.5,
//                        "case recovery" 0.8, "done" 1; each uplift case
//   buckling             entry, "pre-stress" 0, "eigensolve" 0.5, "done" 1
//   modal                entry, "assembly" 0, "eigensolve" / "Lanczos" 0.5
//   second order         every P-Delta iteration ("P-Delta iteration i")
//   corotational         every Newton step; "substep k/n" per load substep
//   arc length           every continuation step ("step s/max_steps")
//   ruin load            every bracketing / bisection probe, and every
//                        Newton step inside a probe (RuinOptions::corot)
//   cable statics        every Newton step; "substep k/n"
//   harmonic, PSD        every sweep frequency ("frequency i/n")
//   influence sweeps     "assembly", "factorization", "adjoint solve", then
//                        every 64 unit-load positions ("ordinates k/n")
//   phased construction  every phase ("phase i/n")
// Thread safety: the token is an atomic flag shared by its copies (cancel()
// from a UI thread, poll from the solver thread); the progress callback runs
// on the solving thread.
#pragma once
#include <functional>
#include <optional>
#include <string>

#include "keta/cancel.h"
#include "keta/core.h"

namespace keta {

// Thrown by a hooked driver whose token fired. A distinct type so callers
// can tell a cancellation from a failure.
class Cancelled : public Error {
 public:
  Cancelled() : Error("analysis cancelled") {}
};

struct SolveHooks {
  std::optional<CancelToken> cancel;                                // none: never cancelled
  std::function<void(double frac, const std::string& phase)> progress;  // empty: no reporting

  // Throw Cancelled if the token fired.
  void check() const {
    if (cancel && cancel->is_cancelled()) throw Cancelled();
  }
  // Report a fraction (clamped to [0, 1]) with a phase label.
  void report(double frac, const std::string& phase) const {
    if (progress) progress(frac < 0 ? 0.0 : frac > 1 ? 1.0 : frac, phase);
  }
};

// Null-safe helpers for the drivers (options carry a nullable pointer).
inline void hooks_check(const SolveHooks* h) {
  if (h) h->check();
}
inline void hooks_report(const SolveHooks* h, double frac, const std::string& phase) {
  if (h) h->report(frac, phase);
}

}  // namespace keta
