// Shared cancellation flag, dependency-free so both the engine hooks
// (keta/hooks.h) and the parametric graph (keta/param/graph.h, SolveCtx) use
// the one type: a token fired from a UI thread cancels the graph walk AND the
// heavy driver running inside a node.
#pragma once
#include <atomic>
#include <memory>

namespace keta {

// Copies observe the same flag.
class CancelToken {
 public:
  CancelToken() : flag_(std::make_shared<std::atomic<bool>>(false)) {}
  void cancel() const { flag_->store(true, std::memory_order_relaxed); }
  bool is_cancelled() const { return flag_->load(std::memory_order_relaxed); }

 private:
  std::shared_ptr<std::atomic<bool>> flag_;
};

}  // namespace keta
