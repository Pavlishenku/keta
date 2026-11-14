#include "keta/parallel.h"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace keta {

namespace {

thread_local bool tl_in_worker = false;

// One job at a time: the calling thread publishes it, the workers and the
// caller take contiguous chunks of indices from an atomic counter.
class Pool {
 public:
  explicit Pool(int n) {
    for (int i = 0; i < n; ++i) workers_.emplace_back([this]() { loop(); });
  }
  ~Pool() {
    {
      std::lock_guard<std::mutex> lk(mu_);
      stop_ = true;
      ++generation_;
    }
    cv_.notify_all();
    for (std::thread& t : workers_) t.join();
  }
  int size() const { return static_cast<int>(workers_.size()); }

  // Every worker takes part in every job (and the caller too): they draw
  // chunks of indices from an atomic counter until none is left; the caller
  // returns once every worker has checked out of the job.
  void run(std::size_t n, const std::function<void(std::size_t)>& body) {
    std::lock_guard<std::mutex> job_lock(job_mu_);  // one job at a time
    errs_.assign(n, nullptr);
    {
      std::lock_guard<std::mutex> lk(mu_);
      body_ = &body;
      n_ = n;
      chunk_ = std::max<std::size_t>(1, n / (8 * static_cast<std::size_t>(size() + 1)));
      next_.store(0);
      pending_ = size();
      ++generation_;
    }
    cv_.notify_all();
    tl_in_worker = true;  // the caller works too; nested calls run inline
    work();
    tl_in_worker = false;
    {
      std::unique_lock<std::mutex> lk(mu_);
      done_cv_.wait(lk, [this]() { return pending_ == 0; });
      body_ = nullptr;
    }
    for (const std::exception_ptr& e : errs_)
      if (e) std::rethrow_exception(e);
  }

 private:
  void work() {
    for (;;) {
      const std::size_t a = next_.fetch_add(chunk_);
      if (a >= n_) return;
      const std::size_t b = std::min(n_, a + chunk_);
      for (std::size_t i = a; i < b; ++i) {
        try {
          (*body_)(i);
        } catch (...) {
          errs_[i] = std::current_exception();
        }
      }
    }
  }
  void loop() {
    tl_in_worker = true;
    std::size_t seen = 0;
    for (;;) {
      {
        std::unique_lock<std::mutex> lk(mu_);
        cv_.wait(lk, [&]() { return generation_ != seen; });
        seen = generation_;
        if (stop_) return;
      }
      work();
      {
        std::lock_guard<std::mutex> lk(mu_);
        if (--pending_ == 0) done_cv_.notify_one();
      }
    }
  }

  std::vector<std::thread> workers_;
  std::mutex mu_, job_mu_;
  std::condition_variable cv_, done_cv_;
  std::size_t generation_ = 0;
  bool stop_ = false;
  const std::function<void(std::size_t)>* body_ = nullptr;
  std::size_t n_ = 0, chunk_ = 1;
  std::atomic<std::size_t> next_{0};
  int pending_ = 0;
  std::vector<std::exception_ptr> errs_;
};

Pool& pool() {
  static Pool p(std::max(1, static_cast<int>(std::thread::hardware_concurrency())) - 1);
  return p;
}

}  // namespace

int resolve_threads(int threads) {
  return threads == 0 ? std::max(1, static_cast<int>(std::thread::hardware_concurrency()))
                      : std::max(1, threads);
}

void parallel_for(std::size_t n, int threads, const std::function<void(std::size_t)>& body) {
  const int want = std::min<int>(resolve_threads(threads), static_cast<int>(std::max<std::size_t>(n, 1)));
  if (want <= 1 || tl_in_worker || pool().size() == 0) {
    std::exception_ptr first;
    for (std::size_t i = 0; i < n; ++i) {
      try {
        body(i);
      } catch (...) {
        if (!first) first = std::current_exception();
      }
    }
    if (first) std::rethrow_exception(first);
    return;
  }
  pool().run(n, body);
}

}  // namespace keta
