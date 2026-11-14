#include "keta/suspension.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace keta {

namespace {

// ---------------------------------------------------------------------------
// Cancellation-free elastic-catenary segment kinematics (see suspension.h (3)).
// ---------------------------------------------------------------------------

// asinh(a)/a, continuous at a = 0 (value 1). Used so the weightless limit of the
// projection integrals is EXACT instead of a 0/0 built from two nearly equal
// asinh's.
double asinhc(double a) {
  if (std::fabs(a) < 1e-8) return 1.0 - a * a / 6.0;  // series, |error| < 1e-33
  return std::asinh(a) / a;
}

struct SegState {
  double H = 0;   // horizontal tension [N] (> 0)
  double V1 = 0;  // vertical tension component at the segment start [N]
  double w = 0;   // weight per unit unstressed length [N/m] (>= 0)
  double EA = 0;  // axial rigidity [N] (> 0)
};

// Horizontal / vertical projection of the first `s` of unstressed length.
void seg_projections(const SegState& g, double s, double& lx, double& lz) {
  double V2 = g.V1 + g.w * s;
  double T1 = std::hypot(g.H, g.V1), T2 = std::hypot(g.H, V2);
  double B = T1 - g.V1 * (g.V1 + V2) / (T1 + T2);  // = H^2/T1 when w = 0
  double A = g.w * s * B / (g.H * g.H);
  lx = (s * B / g.H) * asinhc(A) + g.H * s / g.EA;
  lz = s * (g.V1 + V2) / (T1 + T2) + (s / g.EA) * (g.V1 + g.w * s / 2);
}

// Unstressed length whose horizontal projection is exactly `dx`.
// lx(s) is strictly increasing (dlx/ds = H/T2 + H/EA > 0), so a bracket +
// safeguarded Newton converges unconditionally. Deterministic iteration budget.
double seg_length_for_dx(const SegState& g, double dx) {
  require(dx > 0, "suspension: segment with a non-positive horizontal span");
  auto f = [&](double s) {
    double lx, lz;
    seg_projections(g, s, lx, lz);
    return lx - dx;
  };
  double lo = 0.0, hi = dx;
  for (int i = 0; i < 80 && f(hi) < 0; ++i) hi *= 1.5;
  require(f(hi) >= 0, "suspension: cannot reach the segment span (unstable funicular)");
  double s = 0.5 * (lo + hi);
  for (int it = 0; it < 100; ++it) {
    double r = f(s);
    if (r > 0)
      hi = s;
    else
      lo = s;
    if (std::fabs(r) <= 1e-14 * dx) break;
    double V2 = g.V1 + g.w * s;
    double dfds = g.H / std::hypot(g.H, V2) + g.H / g.EA;
    double sn = (dfds > 0) ? s - r / dfds : 0.5 * (lo + hi);
    s = (sn > lo && sn < hi) ? sn : 0.5 * (lo + hi);
  }
  require(std::isfinite(s) && s > 0, "suspension: segment length solve failed");
  return s;
}

// ---------------------------------------------------------------------------
// The forward march: given the first-segment state (H, V0) walk the funicular
// from the left saddle, applying the vector node jumps (2) + (2').
// ---------------------------------------------------------------------------

struct MarchOut {
  std::vector<double> dz;  // elevation gain per segment [m]
  std::vector<double> l0;  // unstressed length per segment [m]
  std::vector<double> h;   // horizontal tension per segment [N]
  std::vector<double> v1;  // vertical tension component at each segment start [N]
  double z_end = 0;        // elevation of the right saddle relative to the left
  double v_end = 0;        // vertical tension component arriving at the right saddle
};

// `deck_x` / `deck_z` are the hanger deck points RELATIVE to the left saddle
// (x along the plan direction, z elevation difference); empty = every hanger
// vertical, which reproduces the classical constant-H march bit for bit (the
// jump (2') is then never evaluated).
MarchOut march(double H0, double V0, double w, double EA, const std::vector<double>& dx,
               const std::vector<double>& hanger_load, const std::vector<double>& deck_x,
               const std::vector<double>& deck_z) {
  MarchOut m;
  std::size_t ns = dx.size();
  m.dz.resize(ns);
  m.l0.resize(ns);
  m.h.resize(ns);
  m.v1.resize(ns);
  bool inclined = !deck_x.empty();
  double H = H0, V = V0, z = 0, x = 0;
  for (std::size_t k = 0; k < ns; ++k) {
    SegState g{H, V, w, EA};
    m.h[k] = H;
    m.v1[k] = V;
    double s = seg_length_for_dx(g, dx[k]);
    double lx, lz;
    seg_projections(g, s, lx, lz);
    m.l0[k] = s;
    m.dz[k] = lz;
    z += lz;
    x += dx[k];
    V += w * s;  // segment weight
    if (k + 1 < ns) {
      double P = hanger_load[k];
      if (inclined) {
        // Vector jump (2'): the hanger force is aligned with the cable-node ->
        // deck-point axis; its horizontal component P * ddx / (-ddz) transfers
        // horizontal tension across the node.
        double ddx = deck_x[k] - x, ddz = deck_z[k] - z;
        require(ddz < 0, "suspension: hanger " + std::to_string(k) +
                             " would be in compression (its deck point is not below the cable)");
        H -= P * ddx / (-ddz);
        require(H > 0,
                "suspension: unstable funicular (the horizontal tension went non-positive at "
                "hanger " +
                    std::to_string(k) + ")");
      }
      V += P;  // vertical jump (2)
    }
  }
  m.z_end = z;
  m.v_end = V;
  return m;
}

// ---------------------------------------------------------------------------
// Back-stay catenary solves (Lot 3). End 1 of the 2-D frame is the SADDLE,
// local +x pointing from the saddle toward the anchorage.
// ---------------------------------------------------------------------------

// Frictionless saddle: unstressed length such that the back-stay tension AT
// THE SADDLE equals `T_target` (rope-on-pulley equilibrium). Bisection on L0
// in the taut branch — T(saddle) decreases strictly with L0 there — to a
// deterministic budget; every failure is typed.
double backstay_l0_for_saddle_tension(double dx, double dz, double EA, double w,
                                      double T_target) {
  require(T_target > 0, "suspension: back-stay tension target must be positive");
  double chord = std::hypot(dx, dz);
  auto T_of = [&](double l0) {
    Catenary2D c = solve_catenary(dx, dz, l0, EA, w);
    return std::hypot(c.H, c.V1);  // end 1 = saddle
  };
  double l_lo = chord * 0.8, l_hi = chord * 1.2;
  for (int i = 0; i < 200 && T_of(l_lo) < T_target; ++i) l_lo *= 0.98;
  for (int i = 0; i < 200 && T_of(l_hi) > T_target; ++i) l_hi *= 1.02;
  require(T_of(l_lo) >= T_target && T_of(l_hi) <= T_target,
          "suspension: the frictionless-saddle tension equality is unreachable for this "
          "anchorage geometry");
  double l0 = 0.5 * (l_lo + l_hi);
  for (int it = 0; it < 200; ++it) {
    l0 = 0.5 * (l_lo + l_hi);
    double r = T_of(l0) - T_target;
    if (std::fabs(r) <= 1e-12 * T_target) break;
    if (r > 0)
      l_lo = l0;
    else
      l_hi = l0;
  }
  require(std::fabs(T_of(l0) - T_target) <= 1e-9 * T_target,
          "suspension: frictionless-saddle tension equality did not converge");
  return l0;
}

// Solve one back stay. `uh_out` is the plan unit vector from the saddle AWAY
// from the main span; `H_adjacent` the horizontal tension of the adjacent
// main-cable segment; `T_main_saddle` the main-cable tension at the saddle.
BackstayFound solve_backstay(const Vec3& saddle, const BackstaySpec& bs, const Vec3& uh_out,
                             double H_adjacent, double T_main_saddle, double EA, double w) {
  Vec3 rel = bs.anchor - saddle;
  double dxb = rel.x() * uh_out.x() + rel.y() * uh_out.y();
  require(dxb > 1e-9 * (std::fabs(rel.z()) + 1.0),
          "suspension: the back anchorage must lie on the far side of the saddle (away from "
          "the span), not above or inside it");
  double dzb = rel.z();

  BackstayFound bf;
  if (bs.regime == SaddleRegime::Clamped) {
    double Hb = bs.back_H ? *bs.back_H : H_adjacent;
    require(finite(Hb) && Hb > 0, "suspension: back-stay tension target must be positive");
    bf.l0 = catenary_l0_for_tension(dxb, dzb, EA, w, Hb);
  } else {
    require(!bs.back_H,
            "suspension: back_H is a Clamped-regime target — the Frictionless regime "
            "determines the back-stay tension from the tension equality");
    bf.l0 = backstay_l0_for_saddle_tension(dxb, dzb, EA, w, T_main_saddle);
  }
  Catenary2D c = solve_catenary(dxb, dzb, bf.l0, EA, w);
  bf.H = c.H;
  bf.T_saddle = c.T1();
  bf.T_anchor = c.T2();
  Vec3 ez(0, 0, 1);
  bf.saddle_force = c.H * uh_out + c.V1 * ez;         // force ON the saddle node
  bf.anchor_force = -(c.H * uh_out + c.V2() * ez);    // force ON the anchorage
  return bf;
}

}  // namespace

// ---------------------------------------------------------------------------
// Form-finding
// ---------------------------------------------------------------------------

FormFound form_find_cable(const FormFindSpec& spec) {
  require(finite(spec.left) && finite(spec.right), "suspension: non-finite saddle point");
  Vec3 d = spec.right - spec.left;
  double span = std::hypot(d.x(), d.y());
  require(span > 0, "suspension: the two saddles share the same plan position");
  double dz_chord = d.z();
  std::size_t n = spec.hanger_x.size();
  require(n >= 1, "suspension: the main cable needs at least one hanger");
  require(spec.hanger_load.size() == n,
          "suspension: hanger_load and hanger_x must have the same size");
  require(finite(spec.EA) && spec.EA > 0, "suspension: cable EA must be positive");
  require(finite(spec.w) && spec.w >= 0, "suspension: cable weight must be >= 0");
  require(spec.target_sag.has_value() != spec.target_H.has_value(),
          "suspension: give exactly one of target_sag / target_H");
  if (spec.target_sag) require(*spec.target_sag > 0, "suspension: target sag must be positive");
  if (spec.target_H) require(*spec.target_H > 0, "suspension: target H must be positive");
  require(spec.max_iterations > 0 && spec.tol_rel > 0, "suspension: bad solver options");

  double prev = 0;
  double p_total = 0;
  for (std::size_t i = 0; i < n; ++i) {
    require(finite(spec.hanger_x[i]) && spec.hanger_x[i] > prev + 1e-12 &&
                spec.hanger_x[i] < span - 1e-12,
            "suspension: hanger abscissas must be strictly increasing inside the span");
    prev = spec.hanger_x[i];
    require(finite(spec.hanger_load[i]) && spec.hanger_load[i] >= 0,
            "suspension: a hanger load must be finite and >= 0 (a hanger can only pull down)");
    p_total += spec.hanger_load[i];
  }
  require(p_total > 0 || spec.w > 0,
          "suspension: the funicular needs a downward load (all hanger loads and the cable "
          "weight are zero)");

  // Inclined hangers: deck points relative to the left saddle. Empty = vertical.
  std::vector<double> deck_x, deck_z;
  if (!spec.hanger_deck.empty()) {
    require(spec.hanger_deck.size() == n,
            "suspension: hanger_deck must be empty (vertical hangers) or give one deck point "
            "per hanger");
    deck_x.reserve(n);
    deck_z.reserve(n);
    for (const HangerAnchor& hd : spec.hanger_deck) {
      require(finite(hd.x) && finite(hd.z), "suspension: non-finite hanger deck point");
      deck_x.push_back(hd.x);
      deck_z.push_back(hd.z - spec.left.z());
    }
  }

  std::vector<double> dx(n + 1);
  for (std::size_t k = 0; k <= n; ++k) {
    double a = (k == 0) ? 0.0 : spec.hanger_x[k - 1];
    double b = (k == n) ? span : spec.hanger_x[k];
    dx[k] = b - a;
  }
  double x_mid = 0.5 * span;

  // --- initial guess: the parabola of the equivalent UDL --------------------
  double w_total = p_total + spec.w * 1.05 * std::hypot(span, dz_chord);
  double H = spec.target_H ? *spec.target_H : w_total * span / (8.0 * *spec.target_sag);
  require(H > 0 && finite(H), "suspension: degenerate initial horizontal tension");
  double V0 = -0.5 * w_total + H * dz_chord / span;

  // --- residuals ------------------------------------------------------------
  // r1: the march must land on the right saddle.
  // r2: the target (mid-span sag, or the first-segment horizontal tension).
  auto sag_of = [&](const MarchOut& m) {
    // Locate the segment holding the horizontal midpoint and interpolate INSIDE
    // it with the exact catenary (never a chord).
    double x0 = 0, z = 0;
    for (std::size_t k = 0; k <= n; ++k) {
      double x1 = x0 + dx[k];
      if (x_mid <= x1 + 1e-12 || k == n) {
        SegState g{m.h[k], m.v1[k], spec.w, spec.EA};
        double s = seg_length_for_dx(g, std::max(1e-12, x_mid - x0));
        double lx, lz;
        seg_projections(g, s, lx, lz);
        return dz_chord * (x_mid / span) - (z + lz);
      }
      z += m.dz[k];
      x0 = x1;
    }
    return 0.0;
  };
  auto residual = [&](double Hv, double V0v, Eigen::Vector2d& r, MarchOut& m) {
    m = march(Hv, V0v, spec.w, spec.EA, dx, spec.hanger_load, deck_x, deck_z);
    r(0) = (m.z_end - dz_chord) / span;
    if (spec.target_sag)
      r(1) = (sag_of(m) - *spec.target_sag) / span;
    else
      r(1) = (Hv - *spec.target_H) / *spec.target_H;
  };

  Eigen::Vector2d r;
  MarchOut m;
  int it = 0;
  residual(H, V0, r, m);
  for (; it < spec.max_iterations; ++it) {
    if (r.lpNorm<Eigen::Infinity>() <= spec.tol_rel) break;
    // Central-difference Jacobian (the residual itself is exact, so the FD
    // error only affects the convergence RATE, never the converged answer).
    double hH = 1e-7 * H, hV = 1e-7 * std::max(std::fabs(V0), H);
    Eigen::Matrix2d J;
    Eigen::Vector2d rp, rm;
    MarchOut tmp;
    residual(H + hH, V0, rp, tmp);
    residual(H - hH, V0, rm, tmp);
    J.col(0) = (rp - rm) / (2 * hH);
    residual(H, V0 + hV, rp, tmp);
    residual(H, V0 - hV, rm, tmp);
    J.col(1) = (rp - rm) / (2 * hV);
    double det = J.determinant();
    require(std::fabs(det) > 1e-300, "suspension: singular form-finding Jacobian");
    Eigen::Vector2d step = J.fullPivLu().solve(-r);
    require(step.allFinite(), "suspension: non-finite form-finding step");
    // Step limiting: H stays positive and never moves by more than 50 %.
    double scale = 1.0;
    if (std::fabs(step(0)) > 0.5 * H) scale = std::min(scale, 0.5 * H / std::fabs(step(0)));
    H += scale * step(0);
    V0 += scale * step(1);
    require(H > 0, "suspension: unstable funicular (the horizontal tension went non-positive)");
    residual(H, V0, r, m);
  }
  require(r.lpNorm<Eigen::Infinity>() <= std::max(spec.tol_rel * 1e3, 1e-10),
          "suspension: form-finding did not converge (the target sag or tension may be "
          "unreachable for this load set)");

  // --- pack ----------------------------------------------------------------
  FormFound ff;
  ff.H = H;
  ff.EA = spec.EA;
  ff.w = spec.w;
  ff.iterations = it;
  ff.segment_l0 = m.l0;
  ff.segment_h = m.h;
  ff.segment_v1 = m.v1;
  ff.sag = sag_of(m);

  Vec3 uh(d.x() / span, d.y() / span, 0.0);
  Vec3 ez(0, 0, 1);
  ff.nodes.reserve(n + 2);
  ff.nodes.push_back(spec.left);
  double z = 0;
  for (std::size_t k = 0; k < n; ++k) {
    z += m.dz[k];
    ff.nodes.push_back(spec.left + uh * spec.hanger_x[k] + Vec3(0, 0, z));
  }
  ff.nodes.push_back(spec.right);
  // The march's last elevation is r1-converged onto the right saddle; snap the
  // stored node to the EXACT saddle so the assembled model is watertight.

  // Hanger axial tensions and the force each exerts on its cable node. For a
  // vertical hanger the axial force IS the carried load; an inclined one
  // carries P as its vertical component, so its axial force is P |d| / (-dz).
  ff.hanger_tension.resize(n);
  ff.hanger_force.resize(n);
  z = 0;
  for (std::size_t k = 0; k < n; ++k) {
    z += m.dz[k];
    double P = spec.hanger_load[k];
    if (deck_x.empty()) {
      ff.hanger_tension[k] = P;
      ff.hanger_force[k] = -P * ez;
    } else {
      double ddx = deck_x[k] - spec.hanger_x[k], ddz = deck_z[k] - z;
      ff.hanger_tension[k] = P * std::hypot(ddx, ddz) / (-ddz);
      ff.hanger_force[k] = (P / (-ddz)) * (ddx * uh + ddz * ez);
    }
  }

  ff.segment_tmax.resize(n + 1);
  ff.total_l0 = 0;
  for (std::size_t k = 0; k <= n; ++k) {
    double v2 = m.v1[k] + spec.w * m.l0[k];
    ff.segment_tmax[k] = std::max(std::hypot(m.h[k], m.v1[k]), std::hypot(m.h[k], v2));
    ff.total_l0 += m.l0[k];
    require(ff.segment_tmax[k] > 0, "suspension: slack segment (zero tension)");
  }
  ff.saddle_force_left = m.h.front() * uh + Vec3(0, 0, m.v1.front());
  ff.saddle_force_right = -(m.h.back() * uh + Vec3(0, 0, m.v_end));
  return ff;
}

double form_found_elevation(const FormFound& ff, double x) {
  require(ff.nodes.size() >= 2, "suspension: empty form-finding result");
  Vec3 d = ff.nodes.back() - ff.nodes.front();
  double span = std::hypot(d.x(), d.y());
  require(x >= -1e-9 && x <= span + 1e-9, "suspension: elevation query outside the span");
  double x0 = 0, z = ff.nodes.front().z();
  for (std::size_t k = 0; k < ff.segment_l0.size(); ++k) {
    Vec3 dk = ff.nodes[k + 1] - ff.nodes[k];
    double dxk = std::hypot(dk.x(), dk.y());
    if (x <= x0 + dxk + 1e-12 || k + 1 == ff.segment_l0.size()) {
      SegState g{ff.segment_h[k], ff.segment_v1[k], ff.w, ff.EA};
      double t = std::max(0.0, std::min(x - x0, dxk));
      if (t <= 0) return z;
      double s = seg_length_for_dx(g, t);
      double lx, lz;
      seg_projections(g, s, lx, lz);
      return z + lz;
    }
    x0 += dxk;
    z += dk.z();
  }
  return z;
}

// ---------------------------------------------------------------------------
// Model assembly
// ---------------------------------------------------------------------------

double hanger_unstrained_length(double length, double EA, double tension, double w_axial) {
  const std::string who = "hanger_unstrained_length";
  require(finite(length) && length > 0 && finite(EA) && EA > 0, who + ": length and EA finite and > 0");
  require(finite(tension) && tension >= 0 && finite(w_axial), who + ": a finite tension >= 0 and a finite weight");
  // length = L0 (1 + T / EA) + w L0^2 / (2 EA): the root that tends to
  // length / (1 + T / EA) as w -> 0, written without cancellation
  const double a = 1.0 + tension / EA;
  const double disc = a * a + 2.0 * w_axial * length / EA;
  require(disc >= 0, who + ": the tension turns negative along the hanger");
  const double l0 = 2.0 * length / (a + std::sqrt(disc));
  require(tension + w_axial * l0 >= 0, who + ": the tension turns negative along the hanger");
  return l0;
}

SuspensionModel build_suspension(Model& model, const SuspensionSpec& spec) {
  require(spec.cable_material != 0 && spec.cable_section != 0,
          "suspension: the main cable needs a material and a section");
  require(spec.hanger_material != 0 && spec.hanger_section != 0,
          "suspension: the hangers need a material and a section");
  require(spec.hanger_kind != BarKind::Cable,
          "suspension: a hanger cannot be a Cable bar (the funicular models a hanger as a "
          "straight axial link, and the exact catenary rejects a vertical chord) — use Truss "
          "or Frame");
  std::size_t n = spec.form.hanger_x.size();
  require(spec.hanger_deck_nodes.size() == n,
          "suspension: one deck node is needed per hanger");

  // Cable EA / w default to the section+material pair, EXACTLY as the cable
  // driver derives them (cable.cpp resolve_cable_l0), so form-finding and
  // solve_cable_static see the same cable.
  FormFindSpec fs = spec.form;
  const Material& cmat = model.material(spec.cable_material);
  const SectionProps& csec = model.section(spec.cable_section).props;
  if (fs.EA <= 0) fs.EA = cmat.E * csec.A;
  if (fs.w <= 0) fs.w = cmat.rho * csec.A * kGravity;

  // The cable plane: vertical plane through the saddle chord.
  Vec3 d = fs.right - fs.left;
  double span = std::hypot(d.x(), d.y());
  require(span > 0, "suspension: the two saddles share the same plan position");
  Vec3 uh(d.x() / span, d.y() / span, 0.0);
  Vec3 oop(-uh.y(), uh.x(), 0.0);  // out-of-plane unit vector

  // Deck attachment positions are read BEFORE any node/bar is appended (adding
  // entities invalidates the lazy id maps; caching keeps the assembly linear).
  // They drive the hanger DIRECTIONS: `form.hanger_deck` is derived from them,
  // so a deck node vertically under its cable abscissa gives the classical
  // vertical hanger (bit for bit) and an offset one gives an inclined hanger.
  std::vector<Vec3> deck_pos;
  deck_pos.reserve(n);
  for (Id nd : spec.hanger_deck_nodes) deck_pos.push_back(model.node(nd).x);
  fs.hanger_deck.clear();
  fs.hanger_deck.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    Vec3 rel = deck_pos[i] - fs.left;
    double off = rel.x() * oop.x() + rel.y() * oop.y();
    require(std::fabs(off) <= spec.plane_tol * std::max(span, 1.0),
            "suspension: deck node of hanger " + std::to_string(i) +
                " is out of the cable plane — one cable plane per call (build each cable "
                "plane separately)");
    fs.hanger_deck.push_back({rel.x() * uh.x() + rel.y() * uh.y(), deck_pos[i].z()});
  }

  SuspensionModel sm;
  sm.form = form_find_cable(fs);
  sm.deck_nodes = spec.hanger_deck_nodes;
  sm.deck_bars = spec.deck_bars;

  // --- cable nodes ----------------------------------------------------------
  sm.cable_nodes.resize(n + 2, 0);
  auto place = [&](std::size_t k, Id existing) {
    const Vec3& p = sm.form.nodes[k];
    if (existing != 0) {
      const Node& nd = model.node(existing);
      require((nd.x - p).norm() <= spec.saddle_pos_tol * std::max(1.0, p.norm()),
              "suspension: the given saddle node is not at the saddle point of the spec");
      sm.cable_nodes[k] = existing;
      return;
    }
    Id id = model.add_node(p.x(), p.y(), p.z());
    sm.cable_nodes[k] = id;
    if (spec.fix_saddles && (k == 0 || k == n + 1))
      model.add_support(id, {true, true, true, true, true, true});
  };
  place(0, spec.left_node);
  for (std::size_t i = 1; i <= n; ++i) place(i, 0);
  place(n + 1, spec.right_node);

  // --- main-cable segments --------------------------------------------------
  for (std::size_t k = 0; k <= n; ++k) {
    Id b = model.add_bar(sm.cable_nodes[k], sm.cable_nodes[k + 1], spec.cable_material,
                         spec.cable_section, BarKind::Cable);
    model.bars.back().cable = CableProps::with_length(sm.form.segment_l0[k]);
    sm.cable_bars.push_back(b);
  }

  // --- hangers --------------------------------------------------------------
  require(finite(spec.hanger_weight) && spec.hanger_weight >= 0, "suspension: hanger_weight finite and >= 0");
  const double hanger_EA = model.material(spec.hanger_material).E * model.section(spec.hanger_section).props.A;
  require(hanger_EA > 0, "suspension: the hangers need E A > 0");
  for (std::size_t i = 0; i < n; ++i) {
    Id top = sm.cable_nodes[i + 1];
    Id bot = spec.hanger_deck_nodes[i];
    Vec3 pt = sm.form.nodes[i + 1];
    Vec3 pb = deck_pos[i];
    // The march already rejected a compression hanger; keep a defensive check
    // so a hand-built FormFound cannot smuggle one in.
    require(pt.z() - pb.z() > 0, "suspension: hanger " + std::to_string(i) +
                                     " has a non-positive drop (the deck node is not below "
                                     "the cable)");
    sm.hanger_bars.push_back(
        model.add_bar(top, bot, spec.hanger_material, spec.hanger_section, spec.hanger_kind));
    // Unstrained length: the form-found axial tension is the one at the cable
    // node, and falls towards the foot by the axial part of the hanger weight.
    const double len = (pt - pb).norm();
    sm.hanger_l0.push_back(hanger_unstrained_length(len, hanger_EA, sm.form.hanger_tension[i],
                                                    -spec.hanger_weight * (pt.z() - pb.z()) / len));
  }
  // Lack of fit: the bar of geometric length L is fabricated at L0, the
  // imposed strain -(L - L0) / L0 gives it the mean tension EA (L - L0) / L0
  // of the hanger while its ends stay on the form.
  if (spec.hanger_fit_case != 0) {
    LoadCase* fit = nullptr;
    for (LoadCase& c : model.cases)
      if (c.id == spec.hanger_fit_case) fit = &c;
    require(fit != nullptr, "suspension: no load case " + std::to_string(spec.hanger_fit_case) +
                                " to receive the lack of fit of the hangers");
    for (std::size_t i = 0; i < n; ++i) {
      const double len = (sm.form.nodes[i + 1] - deck_pos[i]).norm();
      fit->loads.push_back(Load::prestrain(sm.hanger_bars[i], -(len - sm.hanger_l0[i]) / sm.hanger_l0[i]));
    }
  }

  // --- back stays -----------------------------------------------------------
  // Main-cable tensions at the saddles, for the frictionless equality.
  double T_left = std::hypot(sm.form.segment_h.front(), sm.form.segment_v1.front());
  double v_end = sm.form.segment_v1.back() + sm.form.w * sm.form.segment_l0.back();
  double T_right = std::hypot(sm.form.segment_h.back(), v_end);
  auto add_backstay = [&](const BackstaySpec& bs, std::size_t saddle_ix, const Vec3& uh_out,
                          double H_adjacent, double T_main) {
    require(finite(bs.anchor), "suspension: non-finite back-anchorage point");
    Vec3 rel = bs.anchor - fs.left;
    double off = rel.x() * oop.x() + rel.y() * oop.y();
    require(std::fabs(off) <= spec.plane_tol * std::max(span, 1.0),
            "suspension: the back anchorage is out of the cable plane");
    const Vec3& S = sm.form.nodes[saddle_ix];
    BackstayFound bf = solve_backstay(S, bs, uh_out, H_adjacent, T_main, fs.EA, fs.w);
    Id anchor_id = bs.anchor_node;
    if (anchor_id == 0) {
      anchor_id = model.add_node(bs.anchor.x(), bs.anchor.y(), bs.anchor.z());
      model.add_support(anchor_id, {true, true, true, true, true, true});
    } else {
      require((model.node(anchor_id).x - bs.anchor).norm() <=
                  spec.saddle_pos_tol * std::max(1.0, bs.anchor.norm()),
              "suspension: the given anchorage node is not at the anchor point of the spec");
    }
    Id b = model.add_bar(sm.cable_nodes[saddle_ix], anchor_id, spec.cable_material,
                         spec.cable_section, BarKind::Cable);
    model.bars.back().cable = CableProps::with_length(bf.l0);
    sm.backstay_bars.push_back(b);
    sm.anchor_nodes.push_back(anchor_id);
    return bf;
  };
  if (spec.back_left)
    sm.back_left = add_backstay(*spec.back_left, 0, -uh, sm.form.segment_h.front(), T_left);
  if (spec.back_right)
    sm.back_right =
        add_backstay(*spec.back_right, n + 1, uh, sm.form.segment_h.back(), T_right);

  model.invalidate_index();
  return sm;
}

// ---------------------------------------------------------------------------
// Tower-flexibility coupling
// ---------------------------------------------------------------------------

SuspensionFlexResult build_suspension_on_towers(Model& model, const SuspensionSpec& spec,
                                                Id dead_case,
                                                const SuspensionFlexOptions& opt) {
  require(spec.left_node != 0 && spec.right_node != 0,
          "suspension: the tower coupling needs existing tower-top saddle nodes (left_node "
          "and right_node)");
  require(opt.max_iterations > 0 && opt.saddle_tol > 0,
          "suspension: bad tower-coupling options");
  const Vec3 pL = model.node(spec.left_node).x;
  const Vec3 pR = model.node(spec.right_node).x;
  require((pL - spec.form.left).norm() <= 1e-6 * std::max(1.0, pL.norm()) &&
              (pR - spec.form.right).norm() <= 1e-6 * std::max(1.0, pR.norm()),
          "suspension: form.left / form.right must be the undeformed tower-top positions");
  (void)model.load_case(dead_case);  // typed lookup: the permanent case must exist

  // Fixed point on the saddle displacements (see the header comment for the
  // contraction argument). The trial saddle POINT is the displaced position
  // while the saddle NODE stays at the undeformed one — loosen the position
  // check accordingly, by exactly the displacement magnitude. The hanger
  // cable-attachment abscissas are RE-EXPRESSED in the displaced chord frame
  // so the attachment points stay put in space: a tower lean must not drag
  // the cable clamps sideways off their deck nodes (that would tilt every
  // hanger the same way and feed a spurious horizontal into the deck).
  Vec3 d0 = spec.form.right - spec.form.left;
  double span0 = std::hypot(d0.x(), d0.y());
  require(span0 > 0, "suspension: the two saddles share the same plan position");
  Vec3 uh0(d0.x() / span0, d0.y() / span0, 0.0);
  auto displaced_spec = [&](const Vec3& uL, const Vec3& uR) {
    SuspensionSpec s = spec;
    s.form.left = pL + uL;
    s.form.right = pR + uR;
    Vec3 d1 = s.form.right - s.form.left;
    double span1 = std::hypot(d1.x(), d1.y());
    require(span1 > 0, "suspension: the displaced saddles collapsed in plan");
    Vec3 uh1(d1.x() / span1, d1.y() / span1, 0.0);
    for (double& hx : s.form.hanger_x) {
      Vec3 attach = spec.form.left + uh0 * hx;  // undeformed attachment point
      Vec3 rel = attach - s.form.left;
      hx = rel.x() * uh1.x() + rel.y() * uh1.y();
    }
    if (s.hanger_fit_case == 0) s.hanger_fit_case = dead_case;
    s.saddle_pos_tol = spec.saddle_pos_tol +
                       1.5 * std::max(uL.norm(), uR.norm()) /
                           std::max(1.0, std::min(pL.norm(), pR.norm()));
    return s;
  };

  // The update is ACCELERATED by Aitken / Irons-Tuck dynamic relaxation
  // (Irons & Tuck, IJNME 1969; the standard fixed-point accelerator of
  // partitioned couplings): with back stays the trial model already ties the
  // tower top elastically to ground through the taut back-stay, so the PLAIN
  // iteration contracts only like k_cable/(k_cable + k_tower) (measured ~0.8);
  // the one-parameter secant relaxation restores convergence in a handful of
  // evaluations. Deterministic — theta is clamped and every failure is typed.
  Vec3 uL = Vec3::Zero(), uR = Vec3::Zero();
  Vec6 r_prev = Vec6::Zero();
  double theta = 1.0, gap = 0;
  int it = 0;
  bool converged = false;
  for (; it < opt.max_iterations && !converged; ++it) {
    Model trial = model;
    SuspensionModel tsm = build_suspension(trial, displaced_spec(uL, uR));
    SuspensionResponse resp = suspension_response(trial, tsm, dead_case);
    Vec3 nL =
        resp.raw.result.node_disp[static_cast<std::size_t>(trial.node_index(spec.left_node))]
            .head<3>();
    Vec3 nR =
        resp.raw.result.node_disp[static_cast<std::size_t>(trial.node_index(spec.right_node))]
            .head<3>();
    // Fixed-point residual: solved top position minus the assumed one. Its
    // norm IS the top-equilibrium gap (force gap = k_tower * |r|).
    Vec6 r;
    r << (nL - uL), (nR - uR);
    gap = std::max(r.head<3>().norm(), r.tail<3>().norm());
    if (gap <= opt.saddle_tol) {
      uL = nL;  // finish on the plain update: u = F(u) at convergence
      uR = nR;
      converged = true;
      break;
    }
    if (it > 0) {
      Vec6 dr = r - r_prev;
      double denom = dr.squaredNorm();
      if (denom > 0) theta = -theta * r_prev.dot(dr) / denom;
      theta = std::clamp(theta, 0.05, 20.0);  // keep the secant sane
    }
    r_prev = r;
    uL += theta * r.head<3>();
    uR += theta * r.tail<3>();
  }
  require(converged,
          "suspension: tower coupling did not converge (saddle update " + std::to_string(gap) +
              " m after " + std::to_string(it) + " iterations)");

  SuspensionFlexResult out;
  out.sm = build_suspension(model, displaced_spec(uL, uR));
  out.saddle_disp_left = uL;
  out.saddle_disp_right = uR;
  out.saddle_gap = gap;
  out.iterations = it;
  return out;
}

// ---------------------------------------------------------------------------
// Response
// ---------------------------------------------------------------------------

SuspensionResponse suspension_response(const Model& model, const SuspensionModel& sm, Id case_id,
                                       const CableStaticOptions& opt) {
  SuspensionResponse out;
  out.raw = solve_cable_static(model, case_id, opt);
  const CableStaticResults& r = out.raw;

  out.segment_H.reserve(sm.cable_bars.size());
  for (Id b : sm.cable_bars) {
    std::size_t e = static_cast<std::size_t>(model.bar_index(b));
    out.segment_H.push_back(r.cable_H[e]);
    out.cable_Tmax = std::max(out.cable_Tmax, r.cable_Tmax[e]);
  }
  out.cable_H = out.segment_H.empty() ? 0.0 : out.segment_H.front();

  out.backstay_H.reserve(sm.backstay_bars.size());
  out.backstay_Tmax.reserve(sm.backstay_bars.size());
  for (Id b : sm.backstay_bars) {
    std::size_t e = static_cast<std::size_t>(model.bar_index(b));
    out.backstay_H.push_back(r.cable_H[e]);
    out.backstay_Tmax.push_back(r.cable_Tmax[e]);
  }

  out.hanger_tension.reserve(sm.hanger_bars.size());
  for (Id b : sm.hanger_bars) {
    std::size_t e = static_cast<std::size_t>(model.bar_index(b));
    // Member end forces are local; N = -f(0) is tension-positive (the member
    // pushes back on its own end i when it is stretched).
    out.hanger_tension.push_back(-r.result.bars[e].end_forces_local(0));
  }

  bool first = true;
  for (Id b : sm.deck_bars) {
    std::size_t e = static_cast<std::size_t>(model.bar_index(b));
    for (const Station& st : r.result.bars[e].stations) {
      out.deck_my_max = first ? st.My : std::max(out.deck_my_max, st.My);
      out.deck_my_min = first ? st.My : std::min(out.deck_my_min, st.My);
      first = false;
    }
  }
  first = true;
  for (Id nd : sm.deck_nodes) {
    double dz = r.result.node_disp[static_cast<std::size_t>(model.node_index(nd))](2);
    out.deck_dz_max = first ? dz : std::max(out.deck_dz_max, dz);
    out.deck_dz_min = first ? dz : std::min(out.deck_dz_min, dz);
    first = false;
  }
  return out;
}

SuspensionDelta suspension_live_response(const Model& model, const SuspensionModel& sm,
                                         Id dead_case, Id live_case,
                                         const CableStaticOptions& opt) {
  require(dead_case != live_case, "suspension: dead and live must be two different cases");
  // A cable system is nonlinear: dead + live must be solved as ONE case, never
  // superposed. Build the combined case in a private copy so the caller's model
  // (and every id in `sm`) stays untouched.
  Model combined = model;
  const LoadCase& dead = model.load_case(dead_case);
  const LoadCase& live = model.load_case(live_case);
  Id new_id = 1;
  for (const LoadCase& lc : combined.cases) new_id = std::max(new_id, lc.id + 1);
  LoadCase total;
  total.id = new_id;
  total.name = dead.name + "+" + live.name;
  total.loads = dead.loads;
  total.loads.insert(total.loads.end(), live.loads.begin(), live.loads.end());
  combined.cases.push_back(std::move(total));
  combined.invalidate_index();

  SuspensionDelta out;
  out.dead = suspension_response(model, sm, dead_case, opt);
  out.total = suspension_response(combined, sm, new_id, opt);
  out.hanger_delta.resize(out.dead.hanger_tension.size());
  for (std::size_t i = 0; i < out.hanger_delta.size(); ++i)
    out.hanger_delta[i] = out.total.hanger_tension[i] - out.dead.hanger_tension[i];
  out.cable_H_delta = out.total.cable_H - out.dead.cable_H;
  out.deck_my_max_delta = out.total.deck_my_max - out.dead.deck_my_max;
  out.deck_my_min_delta = out.total.deck_my_min - out.dead.deck_my_min;
  return out;
}

UelmResult solve_uelm(Model& model, const SuspensionSpec& spec, Id dead_case, const UelmOptions& opt) {
  const std::string who = "solve_uelm";
  require(model.load_case(dead_case).id == dead_case, who + ": no dead case");
  UelmResult out;
  Model work = model;
  SuspensionSpec fitted = spec;
  if (fitted.hanger_fit_case == 0) fitted.hanger_fit_case = dead_case;
  out.sm = build_suspension(work, fitted);  // stage 1
  const SuspensionModel& sm = out.sm;
  const std::size_t nn = sm.cable_nodes.size();
  require(nn >= 3, who + ": no hanger node on the main cable");
  // the sag node: the hanger node farthest below the saddle chord
  const Vec3 a = sm.form.nodes.front(), b = sm.form.nodes.back();
  std::size_t sag = 1;
  double best = -1;
  for (std::size_t i = 1; i + 1 < nn; ++i) {
    const Vec3& p = sm.form.nodes[i];
    const double t = (p.x() - a.x()) / (b.x() - a.x());
    const double below = a.z() + t * (b.z() - a.z()) - p.z();
    if (below > best) {
      best = below;
      sag = i;
    }
  }
  // the square target set: z at every hanger node (at its stage-1 place), x at the sag node
  std::vector<TuningTarget> targets;
  for (std::size_t i = 1; i + 1 < nn; ++i) {
    TuningTarget t;
    t.kind = TuningTarget::Kind::Displacement;
    t.id = sm.cable_nodes[i];
    t.comp = 2;
    targets.push_back(t);
  }
  TuningTarget tx;
  tx.kind = TuningTarget::Kind::Displacement;
  tx.id = sm.cable_nodes[sag];
  tx.comp = 0;
  targets.push_back(tx);
  std::vector<TuningUnknown> unknowns;
  for (Id bar : sm.cable_bars) unknowns.push_back(TuningUnknown{bar});
  // stage 2 at the stage-1 lengths: the deviation the correction removes
  {
    const CableStaticResults r = solve_cable_static(work, dead_case, opt.tuning.statics);
    for (const TuningTarget& t : targets)
      out.deviation_stage1 = std::max(out.deviation_stage1,
                                      std::fabs(r.result.node_disp[static_cast<std::size_t>(work.node_index(t.id))](t.comp)));
  }
  const TuningResult tr = tune_cables(work, dead_case, unknowns, targets, opt.tuning);
  out.iterations = tr.iterations;
  out.l0 = tr.x;
  out.state = tr.state;
  for (double v : tr.achieved) out.deviation = std::max(out.deviation, std::fabs(v));
  for (std::size_t j = 0; j < sm.cable_bars.size(); ++j) {
    Bar& bar = work.bars[static_cast<std::size_t>(work.bar_index(sm.cable_bars[j]))];
    CableProps cp = *bar.cable;
    cp.unstressed_length = tr.x[j];
    cp.pretension = 0;
    cp.end_tension = 0;
    bar.cable = cp;
  }
  model = std::move(work);
  return out;
}

CamberResults solve_deck_camber(const Model& model, Id dead_case, const std::vector<std::pair<Id, double>>& targets,
                                const CamberOptions& opt) {
  const std::string who = "solve_deck_camber";
  model.validate();
  require(!targets.empty(), who + ": no target");
  require(opt.max_iterations >= 1 && opt.tol_rel > 0 && finite(opt.tol_rel), who + ": max_iterations >= 1, tol_rel > 0");
  (void)model.load_case(dead_case);
  std::vector<int> ix;
  for (std::size_t i = 0; i < targets.size(); ++i) {
    const auto& [node, z] = targets[i];
    require(finite(z), who + ": the target of node " + std::to_string(node) + " is not finite");
    const int ni = model.node_index(node);
    for (std::size_t k = 0; k < i; ++k)
      require(targets[k].first != node, who + ": node " + std::to_string(node) + " is targeted twice");
    for (const Support& sp : model.supports)
      require(!(sp.node == node && sp.fixed[2]),
              who + ": node " + std::to_string(node) + " is blocked vertically: no camber can move it");
    ix.push_back(ni);
  }
  bool cables = false;
  for (const Bar& b : model.bars) cables |= b.kind == BarKind::Cable;
  CamberResults out;
  {
    Vec3 lo = Vec3::Constant(std::numeric_limits<double>::infinity()), hi = -lo;
    for (const Node& nd : model.nodes) {
      lo = lo.cwiseMin(nd.x);
      hi = hi.cwiseMax(nd.x);
    }
    out.span_scale = (hi - lo).norm();
  }
  const double tol = opt.tol_rel * out.span_scale;
  out.camber.assign(targets.size(), 0.0);
  out.deviation.assign(targets.size(), 0.0);
  for (int pass = 0;; ++pass) {
    Model m = model;
    for (std::size_t i = 0; i < targets.size(); ++i) m.nodes[static_cast<std::size_t>(ix[i])].x.z() += out.camber[i];
    m.invalidate_index();
    CaseResult st;
    if (cables) {
      st = solve_cable_static(m, dead_case, opt.statics).result;
    } else {
      Model m1 = m;
      m1.combos.clear();
      std::erase_if(m1.cases, [&](const LoadCase& c) { return c.id != dead_case; });
      st = solve_static(m1).cases[0];
    }
    out.max_deviation = 0;
    for (std::size_t i = 0; i < targets.size(); ++i) {
      const std::size_t ni = static_cast<std::size_t>(ix[i]);
      out.deviation[i] = m.nodes[ni].x.z() + st.node_disp[ni](2) - targets[i].second;
      out.max_deviation = std::max(out.max_deviation, std::fabs(out.deviation[i]));
    }
    out.trace.push_back(out.max_deviation);
    if (out.max_deviation <= tol) {
      out.iterations = pass;
      out.cambered = std::move(m);
      out.state = std::move(st);
      return out;
    }
    if (pass >= opt.max_iterations) {
      std::string tr;
      for (double d : out.trace) tr += (tr.empty() ? "" : ", ") + std::to_string(d);
      fail(who + ": no convergence after " + std::to_string(pass) + " passes (max deviation trace [" + tr +
           "] m; the hung-deck regime with stiff hangers contracts slowly)");
    }
    for (std::size_t i = 0; i < targets.size(); ++i) out.camber[i] -= out.deviation[i];
  }
}

}  // namespace keta
