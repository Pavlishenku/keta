#include "keta/section_stress.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <set>
#include <thread>
#include <variant>

#include "keta/mesh2d.h"

namespace keta {

namespace {

constexpr int kCircleSegments = 48;  // chords per polygonized circle
constexpr int kFilletSegments = 8;   // chords per quarter-circle root fillet

std::vector<Vec2> rect_ring(double y0, double z0, double y1, double z1) {
  return {Vec2(y0, z0), Vec2(y1, z0), Vec2(y1, z1), Vec2(y0, z1)};
}

// CCW from (+r, 0).
std::vector<Vec2> circle_ring(double r) {
  std::vector<Vec2> out;
  out.reserve(kCircleSegments);
  for (int i = 0; i < kCircleSegments; ++i) {
    const double t = static_cast<double>(i) * 2.0 * kPi / kCircleSegments;
    out.emplace_back(r * std::cos(t), r * std::sin(t));
  }
  return out;
}

double signed_area(const std::vector<Vec2>& r) {
  double a2 = 0;
  for (std::size_t i = 0; i < r.size(); ++i) {
    const Vec2& p = r[i];
    const Vec2& q = r[(i + 1) % r.size()];
    a2 += p.x() * q.y() - q.x() * p.y();
  }
  return a2 / 2;
}

std::vector<Vec2> oriented(std::vector<Vec2> r, bool ccw) {
  if ((signed_area(r) > 0) != ccw) std::reverse(r.begin(), r.end());
  return r;
}

}  // namespace

int target_elements(MeshResolution r) {
  switch (r) {
    case MeshResolution::Coarse: return 1200;
    case MeshResolution::Fine: return 9000;
    case MeshResolution::Standard: break;
  }
  return 3000;
}

std::vector<std::vector<Vec2>> section_rings(const Section& s) {
  struct V {
    std::vector<std::vector<Vec2>> operator()(const ShapeGeneric&) const { return {}; }
    std::vector<std::vector<Vec2>> operator()(const ShapeRect& r) const {
      return {rect_ring(-r.b / 2, -r.h / 2, r.b / 2, r.h / 2)};
    }
    std::vector<std::vector<Vec2>> operator()(const ShapeCircle& c) const { return {circle_ring(c.d / 2)}; }
    std::vector<std::vector<Vec2>> operator()(const ShapeTube& t) const {
      std::vector<Vec2> inner = circle_ring(t.d / 2 - t.t);
      std::reverse(inner.begin(), inner.end());  // hole: CW
      return {circle_ring(t.d / 2), inner};
    }
    std::vector<std::vector<Vec2>> operator()(const ShapeI& i) const {
      // The 12-vertex I outline, plus the four root fillets of a rolled
      // profile (r > 0) as circular arcs: the stresses must use the same
      // section the analysis did (keta's props carry the fillets).
      const double hb = i.b / 2, hh = i.h / 2, hw = i.tw / 2, tf = i.tf;
      const std::vector<Vec2> poly{Vec2(-hb, -hh),     Vec2(hb, -hh),      Vec2(hb, -hh + tf), Vec2(hw, -hh + tf),
                                   Vec2(hw, hh - tf),  Vec2(hb, hh - tf),  Vec2(hb, hh),       Vec2(-hb, hh),
                                   Vec2(-hb, hh - tf), Vec2(-hw, hh - tf), Vec2(-hw, -hh + tf), Vec2(-hb, -hh + tf)};
      if (!(i.r > 0)) return {poly};
      std::vector<Vec2> out;
      const std::size_t n = poly.size();
      for (std::size_t k = 0; k < n; ++k) {
        if (k != 3 && k != 4 && k != 9 && k != 10) {  // the web-flange re-entrant corners
          out.push_back(poly[k]);
          continue;
        }
        const Vec2 c = poly[k];
        const Vec2 a = (c - poly[(k + n - 1) % n]).normalized();  // arriving direction
        const Vec2 b = (poly[(k + 1) % n] - c).normalized();      // leaving direction
        const Vec2 t1 = c - i.r * a, o = t1 + i.r * b;            // 90-degree fillet
        const double th1 = std::atan2(t1.y() - o.y(), t1.x() - o.x());
        double dth = std::atan2((c + i.r * b - o).y(), (c + i.r * b - o).x()) - th1;
        if (dth > kPi) dth -= 2 * kPi;
        if (dth < -kPi) dth += 2 * kPi;
        for (int j = 0; j <= kFilletSegments; ++j) {
          const double th = th1 + dth * j / kFilletSegments;
          out.emplace_back(o.x() + i.r * std::cos(th), o.y() + i.r * std::sin(th));
        }
      }
      return {out};
    }
    std::vector<std::vector<Vec2>> operator()(const ShapeBox& b) const {
      std::vector<Vec2> inner = rect_ring(-(b.b / 2 - b.tw), -(b.h / 2 - b.tf), b.b / 2 - b.tw, b.h / 2 - b.tf);
      std::reverse(inner.begin(), inner.end());
      return {rect_ring(-b.b / 2, -b.h / 2, b.b / 2, b.h / 2), inner};
    }
    std::vector<std::vector<Vec2>> operator()(const ShapeT& t) const {
      // Same centroid as the closed-form T (flange on top, z from the stem foot).
      const double hw = t.h - t.tf;
      const double Af = t.b * t.tf, Aw = hw * t.tw;
      const double zc = (Af * (hw + t.tf / 2) + Aw * hw / 2) / (Af + Aw);
      const double hb = t.b / 2, ht = t.tw / 2;
      return {{Vec2(-ht, -zc), Vec2(ht, -zc), Vec2(ht, hw - zc), Vec2(hb, hw - zc), Vec2(hb, t.h - zc),
               Vec2(-hb, t.h - zc), Vec2(-hb, hw - zc), Vec2(-ht, hw - zc)}};
    }
    std::vector<std::vector<Vec2>> operator()(const ShapePolygon& p) const {
      auto ring = [](const std::vector<Vec3>& v) {
        std::vector<Vec2> r;
        r.reserve(v.size());
        for (const Vec3& q : v) r.emplace_back(q.x(), q.y());
        return r;
      };
      std::vector<std::vector<Vec2>> out{oriented(ring(p.outer), true)};
      for (const auto& h : p.holes) out.push_back(oriented(ring(h), false));
      // Centroid of the holed polygon (holes CW subtract).
      double a = 0, sy = 0, sz = 0;
      for (const auto& r : out)
        for (std::size_t i = 0; i < r.size(); ++i) {
          const Vec2& u = r[i];
          const Vec2& w = r[(i + 1) % r.size()];
          const double cr = u.x() * w.y() - w.x() * u.y();
          a += cr / 2;
          sy += (u.x() + w.x()) * cr / 6;
          sz += (u.y() + w.y()) * cr / 6;
        }
      if (std::fabs(a) > 0)
        for (auto& r : out)
          for (Vec2& q : r) q -= Vec2(sy / a, sz / a);
      return out;
    }
  };
  return std::visit(V{}, s.shape);
}

Vec3 SectionStressEntry::combine(int k, const std::array<double, 6>& a) const {
  Vec3 v = Vec3::Zero();
  for (int i = 0; i < 6; ++i) v += a[static_cast<std::size_t>(i)] * unit[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
  return v;
}

namespace {

// Boundary nodes: the corners and midside of every edge owned by ONE element.
std::vector<int> boundary_nodes(const std::vector<std::array<int, 6>>& tris) {
  std::map<std::pair<int, int>, std::pair<int, int>> edges;  // sorted corners -> (count, mid)
  auto note = [&](int a, int b, int mid) {
    auto& e = edges[{std::min(a, b), std::max(a, b)}];
    ++e.first;
    e.second = mid;
  };
  for (const auto& t : tris) {
    note(t[0], t[1], t[3]);
    note(t[1], t[2], t[4]);
    note(t[2], t[0], t[5]);
  }
  std::vector<int> out;
  for (const auto& [k, v] : edges)
    if (v.first == 1) {
      out.push_back(k.first);
      out.push_back(k.second);
      out.push_back(v.second);
    }
  return out;
}

std::vector<int> pick_candidates(const SectionStressEntry& e) {
  const int nn = static_cast<int>(e.nodes.size());
  std::set<int> set;
  for (int b : boundary_nodes(e.elements)) set.insert(b);
  auto scan = [&](auto field) {
    int imin = 0, imax = 0, itau = 0, ivm = 0;
    double smin = 1e300, smax = -1e300, tmax = -1, vmax = -1;
    for (int k = 0; k < nn; ++k) {
      const Vec3 v = field(k);
      if (v.x() < smin) smin = v.x(), imin = k;
      if (v.x() > smax) smax = v.x(), imax = k;
      const double t2 = v.y() * v.y() + v.z() * v.z();
      if (t2 > tmax) tmax = t2, itau = k;
      const double vm2 = v.x() * v.x() + 3 * t2;
      if (vm2 > vmax) vmax = vm2, ivm = k;
    }
    return std::array<int, 4>{imin, imax, itau, ivm};
  };
  if (nn > 0) {
    for (const auto& f : e.unit) {
      const auto g = scan([&](int k) { return f[static_cast<std::size_t>(k)]; });
      set.insert(g.begin(), g.end());
    }
    // Signed pairs: the |tau| and vm governing node of each a_i +- a_j.
    for (std::size_t i = 0; i < 6; ++i)
      for (std::size_t j = i + 1; j < 6; ++j)
        for (double sgn : {1.0, -1.0}) {
          const auto g = scan([&](int k) {
            return Vec3(e.unit[i][static_cast<std::size_t>(k)] + sgn * e.unit[j][static_cast<std::size_t>(k)]);
          });
          set.insert(g[2]);
          set.insert(g[3]);
        }
  }
  return {set.begin(), set.end()};
}

// The exact solve of a single-material ring set at ~`target` elements. The
// mesh comes from mesh2 (constrained Delaunay + Ruppert refinement, 30 deg),
// 3-5x faster than the exact engine's built-in fallback triangulator for
// the same element count and properties equal to 1e-6 (IPE 300, box, tube,
// rectangle). mesh2's target is an edge LENGTH bound that lands ~2x the
// elements of an equilateral tiling of that edge, hence the sqrt(2). A
// mesh2 failure falls back to the built-in triangulator.
std::shared_ptr<const ExactSection> solve_rings(const std::vector<std::vector<Vec2>>& rings, double nu, int target) {
  ExactInput in;
  in.outline = rings[0];
  in.holes.assign(rings.begin() + 1, rings.end());
  in.materials = {ExactMaterial{"", 1.0, nu}};
  double area = std::fabs(polygon_area(rings[0]));
  for (std::size_t i = 1; i < rings.size(); ++i) area -= std::fabs(polygon_area(rings[i]));
  // Equilateral element of area A / target.
  in.target_h = std::sqrt(4.0 * area / (std::sqrt(3.0) * target));
  Tri6Mesh mesh;
  try {
    mesh2::MeshInput mi;
    mi.outline = in.outline;
    mi.holes = in.holes;
    mi.target_h = in.target_h * std::sqrt(2.0);
    mi.min_angle_deg = 30;
    mesh = tri6_upgrade(mesh2::triangulate(mi), in);
  } catch (const Error&) {
    mesh = mesh_section(in);
  }
  return std::make_shared<const ExactSection>(analyze_exact(mesh, in.materials));
}

std::shared_ptr<const SectionStressEntry> build_entry(std::shared_ptr<const ExactSection> ex) {
  auto e = std::make_shared<SectionStressEntry>();
  e->nodes = ex->mesh.nodes;
  e->elements = ex->mesh.tris;
  const auto fields = ex->nodal_unit_stresses();  // [N, Vy, Vz, My, Mz, T], one sweep
  for (std::size_t i = 0; i < 6; ++i) {
    auto& u = e->unit[i];
    u.reserve(fields[i].size());
    for (const StressState& st : fields[i]) u.emplace_back(st.sig, st.tau_y, st.tau_z);
  }
  e->candidates = pick_candidates(*e);
  // Separable iff sigma ignores Vy, Vz, T and the shear ignores N, My, Mz
  // EXACTLY (true of the Tri6 engine's formulas; checked, never assumed).
  e->separable = true;
  for (std::size_t k = 0; k < e->nodes.size() && e->separable; ++k)
    e->separable = e->unit[1][k].x() == 0 && e->unit[2][k].x() == 0 && e->unit[5][k].x() == 0 &&
                   e->unit[0][k].y() == 0 && e->unit[0][k].z() == 0 && e->unit[3][k].y() == 0 &&
                   e->unit[3][k].z() == 0 && e->unit[4][k].y() == 0 && e->unit[4][k].z() == 0;
  auto& t = e->table;
  for (int k : e->candidates) {
    const std::size_t i = static_cast<std::size_t>(k);
    t.loc.push_back(e->nodes[i]);
    t.sN.push_back(e->unit[0][i].x());
    t.sMy.push_back(e->unit[3][i].x());
    t.sMz.push_back(e->unit[4][i].x());
    t.tyVy.push_back(e->unit[1][i].y());
    t.tyVz.push_back(e->unit[2][i].y());
    t.tyT.push_back(e->unit[5][i].y());
    t.tzVy.push_back(e->unit[1][i].z());
    t.tzVz.push_back(e->unit[2][i].z());
    t.tzT.push_back(e->unit[5][i].z());
  }
  e->solved = std::move(ex);
  return e;
}

}  // namespace

SectionStressEngine SectionStressEngine::build(const Model& model, MeshResolution res) {
  SectionEngineCache cache;
  return build(model, res, cache);
}

SectionStressEngine SectionStressEngine::build(const Model& model, MeshResolution res,
                                               SectionEngineCache& cache, const SolveHooks* hooks) {
  SectionStressEngine eng;
  const int target = target_elements(res);
  std::set<SectionEngineCache::Key> touched;
  std::map<Id, StressFallback> fallback;
  std::set<std::pair<Id, double>> seen;
  // Cache misses are collected, solved in parallel (independent problems,
  // each into its own slot: the result does not depend on the thread
  // count), then published in bar order.
  struct Job {
    SectionEngineCache::Key key;
    std::vector<std::pair<Id, double>> users;  // (section, nu) slots it fills
    Id section = 0;
    std::vector<std::vector<Vec2>> rings;
    double nu = 0;
    std::shared_ptr<const ExactSection> attached;
    std::shared_ptr<const SectionStressEntry> entry;
    std::string error;
  };
  std::vector<Job> jobs;
  std::map<SectionEngineCache::Key, std::size_t> job_of;
  for (const Bar& bar : model.bars) {
    const double nu = model.material(bar.material).nu;
    if (!seen.insert({bar.section, nu}).second) continue;
    hooks_check(hooks);
    const Section& sec = model.section(bar.section);
    std::shared_ptr<const ExactSection> attached = exact_of(sec);
    const bool composite = attached && attached->materials.size() > 1;
    std::vector<std::vector<Vec2>> rings;
    if (!composite) {
      rings = section_rings(sec);
      if (rings.empty()) {
        fallback[sec.id] = {sec.id, StressFallback::Reason::NoGeometry, ""};
        continue;
      }
    }
    SectionEngineCache::Key key;
    if (composite) {
      std::get<3>(key) = attached.get();
    } else {
      for (const auto& r : rings) {
        SectionEngineCache::Ring kr;
        for (const Vec2& p : r) kr.emplace_back(p.x(), p.y());
        std::get<0>(key).push_back(std::move(kr));
      }
      std::get<1>(key) = nu;
      std::get<2>(key) = target;
    }
    touched.insert(key);
    if (auto it = cache.map_.find(key); it != cache.map_.end()) {
      ++cache.stats_.reuses;
      eng.entries_[{bar.section, nu}] = it->second;
      continue;
    }
    if (auto it = job_of.find(key); it != job_of.end()) {  // same content, another section id
      jobs[it->second].users.push_back({bar.section, nu});
      continue;
    }
    job_of[key] = jobs.size();
    Job j;
    j.key = key;
    j.users = {{bar.section, nu}};
    j.section = sec.id;
    j.rings = std::move(rings);
    j.nu = nu;
    j.attached = composite ? attached : nullptr;
    jobs.push_back(std::move(j));
  }

  auto solve_job = [&](Job& j) {
    try {
      j.entry = build_entry(j.attached ? j.attached : solve_rings(j.rings, j.nu, target));
    } catch (const std::exception& e) {
      j.error = e.what();
    }
  };
  std::size_t nt = std::max<std::size_t>(1, std::thread::hardware_concurrency());
  nt = std::min(nt, jobs.size());
  if (nt <= 1) {
    for (Job& j : jobs) {
      hooks_check(hooks);
      solve_job(j);
    }
  } else {
    hooks_check(hooks);
    std::atomic<std::size_t> next{0};
    std::atomic<bool> stop{false};
    std::vector<std::thread> pool;
    for (std::size_t t = 0; t < nt; ++t)
      pool.emplace_back([&] {
        for (std::size_t i; !stop && (i = next++) < jobs.size();) {
          if (hooks && hooks->cancel && hooks->cancel->is_cancelled()) {
            stop = true;
            break;
          }
          solve_job(jobs[i]);
        }
      });
    for (std::thread& th : pool) th.join();
    hooks_check(hooks);  // a cancellation mid-build: nothing published, nothing evicted
  }
  for (Job& j : jobs) {
    if (!j.entry) {
      fallback[j.section] = {j.section, StressFallback::Reason::SolveFailed, j.error};
      continue;
    }
    cache.map_[j.key] = j.entry;
    ++cache.stats_.builds;
    for (const auto& u : j.users) eng.entries_[u] = j.entry;
  }
  // Complete build: keep only the live working set.
  for (auto it = cache.map_.begin(); it != cache.map_.end();)
    it = touched.count(it->first) ? std::next(it) : cache.map_.erase(it);
  for (auto& [id, f] : fallback) eng.fallback_.push_back(std::move(f));
  return eng;
}

const SectionStressEntry* SectionStressEngine::entry(Id section, double nu) const {
  auto it = entries_.find({section, nu});
  return it == entries_.end() ? nullptr : it->second.get();
}

std::vector<Eigen::Vector4d> SectionStressEngine::nodal_field(Id section, double nu,
                                                              const std::array<double, 6>& a) const {
  std::vector<Eigen::Vector4d> out;
  const SectionStressEntry* e = entry(section, nu);
  if (!e) return out;
  out.reserve(e->nodes.size());
  for (int k = 0; k < static_cast<int>(e->nodes.size()); ++k) {
    const Vec3 v = e->combine(k, a);
    out.emplace_back(v.x(), v.y(), v.z(), std::sqrt(v.x() * v.x() + 3 * (v.y() * v.y() + v.z() * v.z())));
  }
  return out;
}

std::optional<StressState> SectionStressEngine::point_stress(Id section, double nu, const std::array<double, 6>& a,
                                                             const Vec2& yz) const {
  const SectionStressEntry* e = entry(section, nu);
  if (!e || e->solved->locate(yz.x(), yz.y()) < 0) return std::nullopt;
  SectionActions sa;
  sa.N = a[0];
  sa.Vy = a[1];
  sa.Vz = a[2];
  sa.My = a[3];
  sa.Mz = a[4];
  sa.T = a[5];
  return e->solved->stress_at(yz.x(), yz.y(), sa);
}

std::array<double, 6> station_actions(const Station& st) { return {st.N, st.Vy, st.Vz, st.My, st.Mz, st.T}; }

namespace {

// The extremes are tracked on SQUARED shear / von Mises (sqrt is monotone, so
// the governing node is the same) and only the winners are square-rooted:
// no sqrt / hypot per candidate.
StationStress exact_station(const Station& st, const SectionStressEntry& e) {
  StationStress out;
  ShearExtreme sh;
  double smin = 1e300, smax = -1e300, t2max = -1, vm2max = -1;
  std::size_t imin = 0, imax = 0, itau = 0, ivm = 0;
  const std::size_t n = e.candidates.size();
  if (e.separable) {
    const auto& t = e.table;
    const double* sN = t.sN.data();
    const double* sMy = t.sMy.data();
    const double* sMz = t.sMz.data();
    const double* yVy = t.tyVy.data();
    const double* yVz = t.tyVz.data();
    const double* yT = t.tyT.data();
    const double* zVy = t.tzVy.data();
    const double* zVz = t.tzVz.data();
    const double* zT = t.tzT.data();
    for (std::size_t k = 0; k < n; ++k) {
      // Same term order as combine(): N, Vy, Vz, My, Mz, T.
      const double sg = st.N * sN[k] + st.My * sMy[k] + st.Mz * sMz[k];
      const double ty = st.Vy * yVy[k] + st.Vz * yVz[k] + st.T * yT[k];
      const double tz = st.Vy * zVy[k] + st.Vz * zVz[k] + st.T * zT[k];
      const double t2 = ty * ty + tz * tz;
      const double vm2 = sg * sg + 3 * t2;
      if (sg < smin) smin = sg, imin = k;
      if (sg > smax) smax = sg, imax = k;
      if (t2 > t2max) t2max = t2, itau = k;
      if (vm2 > vm2max) vm2max = vm2, ivm = k;
    }
  } else {
    const std::array<double, 6> a = station_actions(st);
    for (std::size_t k = 0; k < n; ++k) {
      const Vec3 v = e.combine(e.candidates[k], a);
      const double t2 = v.y() * v.y() + v.z() * v.z();
      const double vm2 = v.x() * v.x() + 3 * t2;
      if (v.x() < smin) smin = v.x(), imin = k;
      if (v.x() > smax) smax = v.x(), imax = k;
      if (t2 > t2max) t2max = t2, itau = k;
      if (vm2 > vm2max) vm2max = vm2, ivm = k;
    }
  }
  if (n == 0) {
    out.shear = sh;
    return out;
  }
  auto loc = [&](std::size_t k) { return e.nodes[static_cast<std::size_t>(e.candidates[k])]; };
  out.sig_min = smin;
  out.sig_max = smax;
  out.sig_min_at = loc(imin);
  out.sig_max_at = loc(imax);
  sh.tau_max = std::sqrt(std::max(t2max, 0.0));
  sh.tau_max_at = loc(itau);
  sh.vm_max = std::sqrt(std::max(vm2max, 0.0));
  sh.vm_max_at = loc(ivm);
  out.shear = sh;
  return out;
}

StationStress fallback_station(const Model& model, const Bar& bar, const Station& st) {
  const Section& sec = model.section(bar.section);
  StationStress out;
  const std::vector<StressAtPoint> pts = station_stresses(model, bar, st);
  if (pts.empty()) {
    out.sig_min = out.sig_max = sec.props.A != 0 ? st.N / sec.props.A : 0.0;
    return out;
  }
  out.sig_min = 1e300;
  out.sig_max = -1e300;
  for (std::size_t i = 0; i < pts.size(); ++i) {
    const Vec2 loc(sec.stress_points[i].y, sec.stress_points[i].z);
    if (pts[i].sigma < out.sig_min) out.sig_min = pts[i].sigma, out.sig_min_at = loc;
    if (pts[i].sigma > out.sig_max) out.sig_max = pts[i].sigma, out.sig_max_at = loc;
  }
  return out;
}

}  // namespace

StressResults recover_stress_results(const Model& model, const StaticResults& results,
                                     const SectionStressEngine& engine, int threads) {
  // Per-bar metadata, resolved once (the model's lazy index is built here, on
  // this thread, so the workers only read).
  std::map<Id, const SectionStressEntry*> entry_of;
  for (const Bar& b : model.bars) entry_of[b.id] = engine.entry(b.section, model.material(b.material).nu);

  // Pre-sized slots: one per (source, bar); sources = cases then combos.
  StressResults out;
  struct Task {
    const BarResult* br;
    BarStress* dst;
  };
  std::vector<Task> tasks;
  std::size_t work = 0;
  auto shape = [&](const std::vector<CaseResult>& src, std::vector<CaseStress>& dst) {
    dst.resize(src.size());
    for (std::size_t c = 0; c < src.size(); ++c) {
      dst[c].id = src[c].id;
      dst[c].name = src[c].name;
      dst[c].bars.resize(src[c].bars.size());
      for (std::size_t b = 0; b < src[c].bars.size(); ++b) {
        tasks.push_back({&src[c].bars[b], &dst[c].bars[b]});
        work += src[c].bars[b].stations.size();
      }
    }
  };
  shape(results.cases, out.cases);
  shape(results.combos, out.combos);

  auto run = [&](std::size_t lo, std::size_t hi) {
    for (std::size_t i = lo; i < hi; ++i) {
      const BarResult& br = *tasks[i].br;
      BarStress& bs = *tasks[i].dst;
      bs.bar = br.bar;
      const SectionStressEntry* e = entry_of.at(br.bar);
      bs.stations.reserve(br.stations.size());
      for (const Station& st : br.stations)
        bs.stations.push_back(e ? exact_station(st, *e) : fallback_station(model, model.bar(br.bar), st));
    }
  };
  int nt = threads > 0 ? threads : static_cast<int>(std::thread::hardware_concurrency());
  // Below ~2000 station-sources (a few ms) a thread costs more than it saves.
  if (work < 2000) nt = 1;
  nt = std::max(1, std::min(nt, static_cast<int>(tasks.size())));
  if (nt == 1) {
    run(0, tasks.size());
  } else {
    std::vector<std::thread> pool;
    std::vector<std::exception_ptr> errors(static_cast<std::size_t>(nt));
    const std::size_t chunk = (tasks.size() + static_cast<std::size_t>(nt) - 1) / static_cast<std::size_t>(nt);
    for (int t = 0; t < nt; ++t) {
      const std::size_t lo = static_cast<std::size_t>(t) * chunk, hi = std::min(tasks.size(), lo + chunk);
      if (lo < hi)
        pool.emplace_back([&, lo, hi, t] {
          try {
            run(lo, hi);
          } catch (...) {
            errors[static_cast<std::size_t>(t)] = std::current_exception();
          }
        });
    }
    for (std::thread& th : pool) th.join();
    for (const std::exception_ptr& e : errors)
      if (e) std::rethrow_exception(e);
  }
  out.fallback_sections = engine.fallback_sections();
  return out;
}

StressEnvelope stress_envelope(const StressResults& stress, const std::string& name) {
  const bool combos = !stress.combos.empty();
  const std::vector<CaseStress>& src = combos ? stress.combos : stress.cases;
  require(!src.empty(), "stress envelope: no sources: at least one case or combination is required");
  const CaseStress& first = src.front();
  for (std::size_t i = 0; i < src.size(); ++i) {
    require(src[i].bars.size() == first.bars.size(),
            "stress envelope: source " + std::to_string(i) + " has " + std::to_string(src[i].bars.size()) +
                " bar(s), source 0 has " + std::to_string(first.bars.size()));
    for (std::size_t b = 0; b < first.bars.size(); ++b) {
      const auto& sts = src[i].bars[b].stations;
      require(sts.size() == first.bars[b].stations.size(),
              "stress envelope: source " + std::to_string(i) + " bar " + std::to_string(first.bars[b].bar) +
                  " has " + std::to_string(sts.size()) + " station(s), source 0 has " +
                  std::to_string(first.bars[b].stations.size()));
      for (std::size_t s = 0; s < sts.size(); ++s) {
        bool bad = !std::isfinite(sts[s].sig_min) || !std::isfinite(sts[s].sig_max);
        if (sts[s].shear) bad |= !std::isfinite(sts[s].shear->tau_max) || !std::isfinite(sts[s].shear->vm_max);
        require(!bad, "stress envelope: source " + std::to_string(i) + " bar " +
                          std::to_string(first.bars[b].bar) + " station " + std::to_string(s) + " is non-finite");
      }
    }
  }
  StressEnvelope env;
  env.name = name;
  env.sources_are_combos = combos;
  for (const CaseStress& c : src) env.source_ids.push_back(c.id);
  for (std::size_t b = 0; b < first.bars.size(); ++b) {
    BarStressEnv be;
    be.bar = first.bars[b].bar;
    for (std::size_t s = 0; s < first.bars[b].stations.size(); ++s) {
      StationStressEnv e;
      e.has_shear = first.bars[b].stations[s].has_shear();
      for (std::size_t i = 0; i < src.size(); ++i) {
        const StationStress& c = src[i].bars[b].stations[s];
        const int k = static_cast<int>(i);
        if (i == 0 || c.sig_min < e.sig_min) e.sig_min = c.sig_min, e.sig_min_src = k;
        if (i == 0 || c.sig_max > e.sig_max) e.sig_max = c.sig_max, e.sig_max_src = k;
        if (e.has_shear) {
          const double t = c.shear ? c.shear->tau_max : 0.0, v = c.shear ? c.shear->vm_max : 0.0;
          if (i == 0 || t > e.tau_max) e.tau_max = t, e.tau_max_src = k;
          if (i == 0 || v > e.vm_max) e.vm_max = v, e.vm_max_src = k;
        }
      }
      be.stations.push_back(e);
    }
    env.bars.push_back(std::move(be));
  }
  return env;
}

}  // namespace keta
