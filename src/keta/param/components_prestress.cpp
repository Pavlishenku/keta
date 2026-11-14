// Tendon, Phase and Phased Analysis components over keta's prestress and
// phasing engines (prestress.h, phasing.h).
//
// Tendon geometry. Profiles are authored over NORMALIZED stations in [0, 1]
// along the chain and scaled to the absolute abscissa s in [0, S], S = the sum
// of the chain line lengths, then stored as keta PolySegments (piecewise cubic
// in s). The eccentricity lives in each bar's own section frame; Assemble
// resolves the chain and flags the bars it runs backwards (Tendon::flips).
//
// Conventions (documented, tested):
//   * profiles are piecewise CUBIC: a Polynomial profile of degree above 3 is
//     a typed reject;
//   * the passage stations must start at 0 and end at 1;
//   * Phased Analysis: a re-jack tensions by the difference with the APPLIED
//     fraction, not with the current effective force; a target below the
//     applied fraction DETENSIONS; a deactivated bar may be re-activated; a
//     K = 1 loss scales a non-constant force profile proportionally.
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>

#include "keta/param/components.h"

namespace keta::param {

namespace {

constexpr double kSmallSlope = 0.4;  // |e'| limit of the small-slope equivalent loads (~23 deg)
constexpr int kSlopeSamples = 32;    // samples per segment (plus the ends)

const Tree& in_at(const std::vector<Tree>& in, std::size_t i) { return input_or_empty(in, i); }

std::vector<double> numbers_of(const Branch<Value>* b) {
  std::vector<double> v;
  if (b)
    for (const auto& o : *b)
      if (o)
        if (auto x = o->as_number()) v.push_back(*x);
  return v;
}
std::vector<Vec3> points_of(const Branch<Value>* b) {
  std::vector<Vec3> v;
  if (b)
    for (const auto& o : *b)
      if (o)
        if (auto p = o->as_point()) v.push_back(*p);
  return v;
}
std::optional<std::int64_t> int_of(const Branch<Value>* b) {
  if (b)
    for (const auto& o : *b)
      if (o)
        if (auto i = o->as_int()) return i;
  return std::nullopt;
}

// The first number of a branch as a finite, strictly positive quantity.
std::optional<double> positive_of(const Branch<Value>* b, const char* what) {
  if (b)
    for (const auto& o : *b)
      if (o)
        if (auto x = o->as_number()) {
          if (std::isfinite(*x) && *x > 0) return *x;
          throw ComponentError(std::string("tendon: ") + what + " = " + Value::number(*x).describe() +
                               " must be finite and strictly positive");
        }
  return std::nullopt;
}

void check_station(double x, const char* what, std::size_t i) {
  if (!std::isfinite(x)) throw ComponentError("tendon: " + std::string(what) + " " + std::to_string(i + 1) + " is not finite");
  if (x < -1e-9 || x > 1.0 + 1e-9)
    throw ComponentError("tendon: " + std::string(what) + " " + std::to_string(i + 1) + " = " +
                         Value::number(x).describe() + " is out of the normalized range [0, 1]");
}

[[noreturn]] void profile_error(const std::string& what) { throw ComponentError("tendon: " + what); }

std::string fmt4(double x) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.4f", x);
  return buf;
}

// Finite coefficients and |e'| <= 0.4, sampled 32 times per segment.
std::vector<PolySegment> validated(std::vector<PolySegment> segs) {
  auto d1 = [](const std::array<double, 4>& c, double x) { return c[1] + x * (2 * c[2] + x * 3 * c[3]); };
  for (const PolySegment& g : segs) {
    bool finite = std::isfinite(g.s0) && std::isfinite(g.s1);
    for (int k = 0; k < 4; ++k) finite = finite && std::isfinite(g.cy[k]) && std::isfinite(g.cz[k]);
    if (!finite) profile_error("tendon profile received a non-finite value");
    const double span = g.s1 - g.s0;
    for (int i = 0; i <= kSlopeSamples; ++i) {
      const double x = span * i / kSlopeSamples;
      const double slope = std::hypot(d1(g.cy, x), d1(g.cz, x));
      if (!std::isfinite(slope)) profile_error("tendon profile received a non-finite value");
      if (slope > kSmallSlope)
        profile_error("tendon profile slope " + fmt4(slope) + " exceeds the small-slope limit " + fmt4(kSmallSlope));
    }
  }
  return segs;
}

// (s, e_y, e_z) vertices: strictly increasing, finite, spanning [0, S].
void check_vertices(const std::vector<std::array<double, 3>>& pts, double S) {
  for (const auto& p : pts)
    if (!std::isfinite(p[0]) || !std::isfinite(p[1]) || !std::isfinite(p[2]))
      profile_error("tendon profile received a non-finite value");
  for (std::size_t i = 1; i < pts.size(); ++i)
    if (pts[i][0] <= pts[i - 1][0])
      profile_error("tendon profile abscissae must be strictly increasing (at index " + std::to_string(i) + ")");
  if (std::fabs(pts.front()[0]) > 1e-9 * S || std::fabs(pts.back()[0] - S) > 1e-9 * S)
    profile_error("the passage stations must run from 0 to 1 (keta profiles span the whole chain)");
}

std::vector<PolySegment> polyline(const std::vector<std::array<double, 3>>& pts, double S) {
  check_vertices(pts, S);
  std::vector<PolySegment> segs;
  for (std::size_t i = 0; i + 1 < pts.size(); ++i) {
    const double h = pts[i + 1][0] - pts[i][0];
    PolySegment g;
    g.s0 = pts[i][0];
    g.s1 = pts[i + 1][0];
    g.cy = {pts[i][1], (pts[i + 1][1] - pts[i][1]) / h, 0, 0};
    g.cz = {pts[i][2], (pts[i + 1][2] - pts[i][2]) / h, 0, 0};
    segs.push_back(g);
  }
  segs.front().s0 = 0;
  segs.back().s1 = S;
  return validated(std::move(segs));
}

// Natural cubic spline second derivatives (Thomas algorithm).
std::vector<double> spline_moments(const std::vector<double>& s, const std::vector<double>& v) {
  const std::size_t n = s.size();
  std::vector<double> m(n, 0.0);
  if (n < 3) return m;
  std::vector<double> h(n - 1);
  for (std::size_t i = 0; i + 1 < n; ++i) h[i] = s[i + 1] - s[i];
  const std::size_t ni = n - 2;
  std::vector<double> a(ni), b(ni), c(ni), d(ni);
  for (std::size_t k = 0; k < ni; ++k) {
    const std::size_t i = k + 1;
    a[k] = h[i - 1];
    b[k] = 2 * (h[i - 1] + h[i]);
    c[k] = h[i];
    d[k] = 6 * ((v[i + 1] - v[i]) / h[i] - (v[i] - v[i - 1]) / h[i - 1]);
  }
  for (std::size_t k = 1; k < ni; ++k) {
    const double w = a[k] / b[k - 1];
    b[k] -= w * c[k - 1];
    d[k] -= w * d[k - 1];
  }
  std::vector<double> x(ni);
  x[ni - 1] = d[ni - 1] / b[ni - 1];
  for (std::size_t k = ni - 1; k-- > 0;) x[k] = (d[k] - c[k] * x[k + 1]) / b[k];
  for (std::size_t k = 0; k < ni; ++k) m[k + 1] = x[k];
  return m;
}

std::vector<PolySegment> spline(const std::vector<std::array<double, 3>>& pts, double S) {
  check_vertices(pts, S);
  std::vector<double> ss, ys, zs;
  for (const auto& p : pts) {
    ss.push_back(p[0]);
    ys.push_back(p[1]);
    zs.push_back(p[2]);
  }
  const auto my = spline_moments(ss, ys), mz = spline_moments(ss, zs);
  auto coeffs = [](double vi, double vj, double mi, double mj, double h) {
    return std::array<double, 4>{vi, (vj - vi) / h - mi * h / 3 - mj * h / 6, mi / 2, (mj - mi) / (6 * h)};
  };
  std::vector<PolySegment> segs;
  for (std::size_t i = 0; i + 1 < ss.size(); ++i) {
    const double h = ss[i + 1] - ss[i];
    PolySegment g;
    g.s0 = ss[i];
    g.s1 = ss[i + 1];
    g.cy = coeffs(ys[i], ys[i + 1], my[i], my[i + 1], h);
    g.cz = coeffs(zs[i], zs[i + 1], mz[i], mz[i + 1], h);
    segs.push_back(g);
  }
  segs.front().s0 = 0;
  segs.back().s1 = S;
  return validated(std::move(segs));
}

// e(s) = e0 + lin (s/S) + quad (s/S)^2 through s = 0, S/2, S.
std::vector<PolySegment> parabola(const std::vector<Vec3>& p, double S) {
  for (const Vec3& q : p)
    if (!std::isfinite(q.y()) || !std::isfinite(q.z())) profile_error("tendon profile received a non-finite value");
  auto comp = [S](double e0, double em, double e1) {
    const double lin = -3 * e0 + 4 * em - e1, quad = 2 * (e0 - 2 * em + e1);
    return std::array<double, 4>{e0, lin / S, quad / (S * S), 0.0};
  };
  PolySegment g;
  g.s0 = 0;
  g.s1 = S;
  g.cy = comp(p[0].y(), p[1].y(), p[2].y());
  g.cz = comp(p[0].z(), p[1].z(), p[2].z());
  return validated({g});
}

// Free polynomial in xi = s / S, ascending powers, converted to the s basis.
std::vector<PolySegment> polynomial(const std::vector<double>& cy, const std::vector<double>& cz, double S) {
  const std::size_t degree = std::max(cy.size(), cz.size()) - (std::max(cy.size(), cz.size()) > 0 ? 1 : 0);
  if (degree > 8) profile_error("tendon polynomial degree " + std::to_string(degree) + " exceeds the maximum 8");
  if (degree > 3)
    profile_error("tendon polynomial degree " + std::to_string(degree) +
                  " exceeds 3: keta's tendon profiles are piecewise cubic");
  PolySegment g;
  g.s0 = 0;
  g.s1 = S;
  for (std::size_t k = 0; k < 4; ++k) {
    g.cy[k] = k < cy.size() ? cy[k] / std::pow(S, static_cast<double>(k)) : 0.0;
    g.cz[k] = k < cz.size() ? cz[k] / std::pow(S, static_cast<double>(k)) : 0.0;
  }
  return validated({g});
}

// Absolute (s, P) breakpoints from the force values and normalized stations.
std::vector<std::pair<double, double>> force_profile(const std::vector<double>& f, const std::vector<double>& st,
                                                     double S) {
  if (f.empty()) throw ComponentError("tendon: no force values (wire F)");
  std::vector<std::pair<double, double>> pts;
  if (f.size() == 1) {
    pts.emplace_back(0.0, f[0]);
  } else if (st.empty()) {
    for (std::size_t i = 0; i < f.size(); ++i)
      pts.emplace_back(S * static_cast<double>(i) / static_cast<double>(f.size() - 1), f[i]);
  } else {
    if (st.size() != f.size())
      throw ComponentError("tendon: " + std::to_string(f.size()) + " force value(s) but " + std::to_string(st.size()) +
                           " station(s) (F and S must be equal length)");
    for (std::size_t i = 0; i < st.size(); ++i) {
      check_station(st[i], "force station", i);
      pts.emplace_back(st[i] * S, f[i]);
    }
  }
  for (std::size_t i = 0; i < pts.size(); ++i) {
    if (!std::isfinite(pts[i].first) || !std::isfinite(pts[i].second) || pts[i].second < 0)
      profile_error("tendon profile received a non-finite value");
    if (i > 0 && pts[i].first <= pts[i - 1].first)
      profile_error("tendon profile abscissae must be strictly increasing (at index " + std::to_string(i) + ")");
  }
  return pts;
}

// Chain mean of a piecewise-linear profile over [0, S], ends held constant.
double mean_over(const std::vector<std::pair<double, double>>& pts, double S) {
  ForceProfile fp;
  for (const auto& [s, p] : pts) {
    fp.s.push_back(s);
    fp.P.push_back(p);
  }
  if (S <= 0) return pts.front().second;
  std::vector<double> xs{0.0, S};
  for (const auto& pt : pts)
    if (pt.first > 0 && pt.first < S) xs.push_back(pt.first);
  std::sort(xs.begin(), xs.end());
  double integral = 0;
  for (std::size_t i = 0; i + 1 < xs.size(); ++i) integral += 0.5 * (fp.at(xs[i]) + fp.at(xs[i + 1])) * (xs[i + 1] - xs[i]);
  return integral / S;
}

std::string canonical_mode(const std::string& m) {
  std::string k;
  for (char c : m)
    if (!std::isspace(static_cast<unsigned char>(c))) k.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
  if (k == "polyline") return "Polyline";
  if (k == "parabola") return "Parabola";
  if (k == "spline") return "Spline";
  if (k == "polynomial") return "Polynomial";
  return {};
}
bool is_pre(const std::string& k) {
  std::string x;
  for (char c : k)
    if (!std::isspace(static_cast<unsigned char>(c))) x.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
  return x == "pre" || x == "pretension" || x == "pre-tension";
}

std::vector<Line> lines_of(const Tree& t) {
  std::vector<Line> out;
  for (const auto& kv : t.branches())
    for (const auto& o : kv.second)
      if (o) {
        if (const Line* l = o->as_line()) out.push_back(*l);
        else if (const auto* b = o->get<BeamSpec>()) out.push_back(b->line);
      }
  return out;
}

}  // namespace

Tendon tendon_record(const TendonSpec& spec) {
  Tendon t;
  t.name = spec.name.value_or("");
  t.profile = spec.profile;
  t.force_profile = spec.force;
  double S = spec.profile.empty() ? 0.0 : spec.profile.back().s1;
  t.jack_force = mean_over(spec.force, S);
  t.losses.mu = t.losses.k = t.losses.slip = 0;  // all inside the authored force
  if (spec.ep && spec.ap) {
    t.losses.Ep = *spec.ep;
    t.losses.Ap = *spec.ap;
  } else {
    t.losses.Ap = 0;  // lossless: the losses need both E_p and A_p
  }
  t.case_id = spec.case_id;
  return t;
}

// ============================ Tendon ============================

ComponentInfo TendonNode::info() const {
  return {"keta.tendon", "Tendon", "Model", "Prestressing tendon: bar chain + eccentricity + force profile"};
}
std::vector<PinSpec> TendonNode::inputs() const {
  return {PinSpec::list("L", "beam", "Chain bars/lines (each a bar; authored travel order)"),
          PinSpec::list("P", "point", "Profile passage points (x = normalized station [0,1], y = e_y m, z = e_z m)").opt(),
          PinSpec::list("F", "number", "Tendon force values (N, P >= 0)"),
          PinSpec::list("S", "number", "Force stations (normalized [0,1]; empty = equally spaced; single F = constant)").opt(),
          PinSpec::list("Y", "number", "e_y polynomial coeffs (xi-basis; Polynomial mode)").opt(),
          PinSpec::list("Z", "number", "e_z polynomial coeffs (xi-basis; Polynomial mode)").opt(),
          PinSpec::item("C", "int", "Load-case id (default 1)").opt(),
          PinSpec::item("EP", "number", "Tendon modulus E_p (Pa); staged-tensioning losses").opt(),
          PinSpec::item("AP", "number", "Tendon area A_p (m^2); staged-tensioning losses").opt()};
}
std::vector<PinSpec> TendonNode::outputs() const { return {PinSpec::item("T", "tendon", "The tendon")}; }
Outputs<Value> TendonNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const std::string m = canonical_mode(mode);
  if (m.empty())
    throw ComponentError("tendon: unknown profile mode \"" + mode + "\" (expected Polyline, Parabola, Spline or Polynomial)");
  std::optional<std::string> nm;
  {
    std::string t = name;
    t.erase(0, t.find_first_not_of(" \t\n\r"));
    t.erase(t.find_last_not_of(" \t\n\r") + 1);
    if (!t.empty()) nm = t;
  }
  if (!has_values(in_at(in, 0))) throw ComponentError::missing("L");

  std::vector<BranchIndex> ix;
  for (std::size_t p = 1; p <= 8; ++p) ix.emplace_back(in_at(in, p));
  auto at = [&](std::size_t pin, std::size_t bi) { return ix[pin - 1].at(bi); };

  Tree out;
  std::size_t bi = 0;
  for (const auto& [path, slots] : in_at(in, 0).branches()) {
    TendonSpec t;
    for (const auto& o : slots)
      if (o) {
        if (const Line* l = o->as_line()) t.lines.push_back(*l);
        else if (const auto* b = o->get<BeamSpec>()) t.lines.push_back(b->line);
      }
    if (t.lines.empty()) throw ComponentError("tendon: a chain branch has no lines (wire the bar chain into L)");
    double S = 0;
    for (std::size_t i = 0; i < t.lines.size(); ++i) {
      const double len = t.lines[i].length();
      if (!(std::isfinite(len) && len > 0))
        throw ComponentError("tendon: chain line " + std::to_string(i + 1) + " has zero or non-finite length");
      S += len;
    }

    const std::vector<Vec3> pts = points_of(at(1, bi));
    if (m == "Polynomial") {
      t.profile = polynomial(numbers_of(at(4, bi)), numbers_of(at(5, bi)), S);
    } else if (m == "Parabola") {
      if (pts.size() != 3)
        throw ComponentError("tendon: Parabola mode needs exactly 3 passage points (e at s=0, S/2, S), got " +
                             std::to_string(pts.size()));
      for (std::size_t i = 0; i < 3; ++i) check_station(pts[i].x(), "passage station", i);
      t.profile = parabola(pts, S);
    } else {
      const std::size_t need = m == "Polyline" ? 2 : 3;
      if (pts.size() < need)
        throw ComponentError("tendon: " + m + " mode needs at least " + std::to_string(need) +
                             " passage points, got " + std::to_string(pts.size()));
      std::vector<std::array<double, 3>> v;
      for (std::size_t i = 0; i < pts.size(); ++i) {
        check_station(pts[i].x(), "passage station", i);
        v.push_back({pts[i].x() * S, pts[i].y(), pts[i].z()});
      }
      t.profile = m == "Polyline" ? polyline(v, S) : spline(v, S);
    }
    t.force = force_profile(numbers_of(at(2, bi)), numbers_of(at(3, bi)), S);
    t.case_id = static_cast<Id>(int_of(at(6, bi)).value_or(1));
    t.ep = positive_of(at(7, bi), "tendon modulus E_p");
    t.ap = positive_of(at(8, bi), "tendon area A_p");
    t.pretension = is_pre(kind);
    t.name = nm;
    out.append(path, Value::of(std::move(t)));
    ++bi;
  }
  return {out};
}
json::Value TendonNode::params() const {
  json::Value p = json::Value::object();
  p.set("mode", mode);
  p.set("kind", kind);
  p.set("name", name);
  return p;
}
void TendonNode::set_params(const json::Value& p) {
  // Normalized to the canonical spelling; an unknown mode is ignored.
  if (const auto* x = p.get("mode"); x && x->is_string())
    if (auto c = canonical_mode(x->as_string()); !c.empty()) mode = c;
  if (const auto* x = p.get("kind"); x && x->is_string()) kind = is_pre(x->as_string()) ? "Pre" : "Post";
  if (const auto* x = p.get("name"); x && x->is_string()) name = x->as_string();
}

// ============================ Phase ============================

ComponentInfo PhaseNode::info() const {
  return {"keta.phase", "Phase", "Model", "One construction stage: activate/remove bars & supports, apply cases, jack tendons"};
}
std::vector<PinSpec> PhaseNode::inputs() const {
  return {PinSpec::item("N", "text", "Phase name").opt(),
          PinSpec::list("A", "beam", "Bars activated this phase").opt(),
          PinSpec::list("D", "beam", "Bars removed this phase").opt(),
          PinSpec::list("S", "support", "Supports added this phase").opt(),
          PinSpec::list("X", "support", "Supports removed this phase").opt(),
          PinSpec::list("C", "int", "Load cases applied this phase").opt(),
          PinSpec::list("T", "tendon", "Tendons jacked (simultaneous)").opt(),
          PinSpec::list("F", "number", "Jack fraction per tendon (repeat-last, default 1.0)").opt()};
}
std::vector<PinSpec> PhaseNode::outputs() const { return {PinSpec::item("PH", "phase", "The phase")}; }
Outputs<Value> PhaseNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  PhaseSpec ph;
  if (auto t = first_text(in_at(in, 0))) {
    std::string x = *t;
    x.erase(0, x.find_first_not_of(" \t\n\r"));
    x.erase(x.find_last_not_of(" \t\n\r") + 1);
    if (!x.empty()) ph.name = x;
  }
  ph.activate = lines_of(in_at(in, 1));
  ph.deactivate = lines_of(in_at(in, 2));
  for (const Value& v : in_at(in, 3).non_nulls())
    if (const auto* s = v.get<SupportSpec>()) ph.add_supports.push_back(*s);
  for (const Value& v : in_at(in, 4).non_nulls())
    if (const auto* s = v.get<SupportSpec>()) ph.remove_supports.push_back(*s);
  for (const Value& v : in_at(in, 5).non_nulls())
    if (auto i = v.as_int()) ph.cases.push_back(*i);
  std::vector<double> fracs;
  for (const Value& v : in_at(in, 7).non_nulls())
    if (auto x = v.as_number()) fracs.push_back(*x);
  std::size_t i = 0;
  for (const Value& v : in_at(in, 6).non_nulls()) {
    const auto* t = v.get<TendonSpec>();
    if (!t) continue;
    const double f = fracs.empty() ? 1.0 : (i < fracs.size() ? fracs[i] : fracs.back());
    if (!std::isfinite(f) || f < 0)
      throw ComponentError("phase: jack fraction " + std::to_string(i + 1) + " = " + Value::number(f).describe() +
                           " must be finite and >= 0");
    ph.tension.emplace_back(*t, f);
    ++i;
  }
  return {single(Value::of(std::move(ph)))};
}

// ============================ Phased Analysis ============================

namespace {

std::string at_point(const Vec3& p) {
  char buf[96];
  std::snprintf(buf, sizeof buf, "(%.3f, %.3f, %.3f)", p.x(), p.y(), p.z());
  return buf;
}
std::string tol_str(double tol) { return Value::number(tol).describe(); }

// The nearest node within tol carrying a model support with the same blocked
// pattern (a typed reject when none matches).
Id support_node(const Model& m, const SupportSpec& s, double tol, std::size_t phase_no) {
  const keta::Node* best = nullptr;
  const keta::Node* nearest = nullptr;
  for (const keta::Node& n : m.nodes) {
    const double d = (n.x - s.point).norm();
    if (d > tol) continue;
    if (!nearest || d < (nearest->x - s.point).norm()) nearest = &n;
    const bool match = std::any_of(m.supports.begin(), m.supports.end(),
                                   [&](const Support& x) { return x.node == n.id && x.fixed == s.blocked; });
    if (match && (!best || d < (best->x - s.point).norm())) best = &n;
  }
  if (best) return best->id;
  if (!nearest)
    throw ComponentError("phase " + std::to_string(phase_no) + ": support at " + at_point(s.point) +
                         " lands on no node (tolerance " + tol_str(tol) + ")");
  throw ComponentError("phase " + std::to_string(phase_no) + ": support at " + at_point(s.point) + " matches node " +
                       std::to_string(nearest->id) + " but no model support there has the same blocked pattern");
}

// Structural tendon match: the record Assemble built from an equal spec.
bool tendon_matches(const Tendon& t, const Tendon& want) {
  if (t.name != want.name || t.case_id != want.case_id || t.force_profile != want.force_profile ||
      t.losses.Ep != want.losses.Ep || t.losses.Ap != want.losses.Ap || t.profile.size() != want.profile.size())
    return false;
  for (std::size_t i = 0; i < t.profile.size(); ++i) {
    const PolySegment &a = t.profile[i], &b = want.profile[i];
    if (a.s0 != b.s0 || a.s1 != b.s1 || a.cy != b.cy || a.cz != b.cz) return false;
    if (a.exact != b.exact || a.R != b.R || a.th != b.th || !(a.C == b.C) || !(a.U == b.U) ||
        !(a.D == b.D))
      return false;
    for (std::size_t k = 0; k < 4; ++k)
      if (!(a.p[k] == b.p[k])) return false;
  }
  return true;
}

}  // namespace

ComponentInfo PhasedAnalyzeNode::info() const {
  return {"keta.phased", "Phased Analysis", "Analyze", "Staged construction-phasing analysis (incremental superposition)"};
}
std::vector<PinSpec> PhasedAnalyzeNode::inputs() const {
  return {PinSpec::item("M", "model", "The model"), PinSpec::list("P", "phase", "The phases, in construction order")};
}
std::vector<PinSpec> PhasedAnalyzeNode::outputs() const {
  return {PinSpec::item("R", "results", "Cumulated results (one pseudo-case per phase)"),
          PinSpec::item("I", "results", "Incremental results (one pseudo-case per phase)"),
          PinSpec::item("P", "phased", "The phased analysis (per-tendon forces)")};
}
Outputs<Value> PhasedAnalyzeNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  ModelRef model;
  for (const Value& v : in_at(in, 0).non_nulls())
    if (const auto* m = v.get<ModelRef>(); m && *m) {
      model = *m;
      break;
    }
  if (!model) throw ComponentError::missing("M");
  std::vector<PhaseSpec> specs;
  for (const Value& v : in_at(in, 1).non_nulls())
    if (const auto* p = v.get<PhaseSpec>()) specs.push_back(*p);
  if (specs.empty()) throw ComponentError::missing("P");
  const double tol = tolerance;

  // ---- resolution: geometry / specs -> ids (typed rejects) ----
  std::set<Id> consumed;
  std::vector<std::pair<const TendonSpec*, Id>> assigned;
  PhasePlan plan;
  std::set<Id> activated_support_nodes;
  std::vector<std::vector<Id>> add_nodes(specs.size());
  for (std::size_t pi = 0; pi < specs.size(); ++pi) {
    const PhaseSpec& ps = specs[pi];
    const std::size_t no = pi + 1;
    Phase ph;
    ph.name = ps.name.value_or("");
    auto bars_of = [&](const std::vector<Line>& lines, const char* what) {
      std::vector<Id> ids;
      for (std::size_t i = 0; i < lines.size(); ++i) {
        auto hit = find_bar(*model, lines[i], tol);
        if (!hit)
          throw ComponentError("phase " + std::to_string(no) + ": " + what + " bar line " + std::to_string(i + 1) +
                               " matches no bar in the model (tolerance " + tol_str(tol) + ")");
        ids.push_back(model->bars[hit->index].id);
      }
      return ids;
    };
    ph.activate_bars = bars_of(ps.activate, "activated");
    ph.deactivate_bars = bars_of(ps.deactivate, "removed");
    for (const SupportSpec& s : ps.add_supports) {
      Id n = support_node(*model, s, tol, no);
      add_nodes[pi].push_back(n);
      activated_support_nodes.insert(n);
    }
    for (const SupportSpec& s : ps.remove_supports) ph.remove_support_nodes.push_back(support_node(*model, s, tol, no));
    ph.cases = ps.cases;
    for (const auto& [spec, frac] : ps.tension) {
      // Re-jack: the identical spec (same wire) targets the same record.
      std::optional<Id> id;
      for (auto it = assigned.rbegin(); it != assigned.rend(); ++it)
        if (*it->first == spec) {
          id = it->second;
          break;
        }
      if (!id) {
        const Tendon want = tendon_record(spec);
        for (const Tendon& t : model->tendons)
          if (!consumed.count(t.id) && tendon_matches(t, want)) {
            consumed.insert(t.id);
            assigned.emplace_back(&spec, t.id);
            id = t.id;
            break;
          }
      }
      if (!id)
        throw ComponentError("phase " + std::to_string(no) + ": tendon" + (spec.name ? " \"" + *spec.name + "\"" : "") +
                             " not found in model — wire the same tendon into Assemble T and Phase T");
      ph.tension.push_back({*id, frac});
    }
    plan.phases.push_back(std::move(ph));
  }

  // Supports named in an activation list are inactive until then: the working
  // model starts without them and each phase adds back the node's supports.
  Model work = *model;
  work.supports.clear();
  for (const Support& s : model->supports)
    if (!activated_support_nodes.count(s.node)) work.supports.push_back(s);
  std::unordered_map<Id, std::vector<const Support*>> supports_of;  // per node, in order
  for (const Support& s : model->supports) supports_of[s.node].push_back(&s);
  for (std::size_t pi = 0; pi < specs.size(); ++pi)
    for (Id n : add_nodes[pi])
      for (const Support* s : supports_of[n]) plan.phases[pi].add_supports.push_back(*s);
  work.invalidate_index();

  // ---- plan validation (every defect at once) ----
  {
    std::vector<std::string> errs;
    std::set<Id> ever, active_bars, active_sup;
    for (const Phase& ph : plan.phases)
      ever.insert(ph.activate_bars.begin(), ph.activate_bars.end());
    for (const Bar& b : model->bars)
      if (!ever.count(b.id)) active_bars.insert(b.id);
    for (const Support& s : work.supports) active_sup.insert(s.node);
    std::map<Id, Id> tendon_case;
    for (const Tendon& t : model->tendons)
      if (t.case_id != 0) tendon_case[t.case_id] = t.id;
    std::set<Id> live;
    for (std::size_t pi = 0; pi < plan.phases.size(); ++pi) {
      const Phase& ph = plan.phases[pi];
      const std::string P = "phase " + std::to_string(pi + 1) + ": ";
      for (Id b : ph.activate_bars)
        if (!active_bars.insert(b).second)
          errs.push_back(P + "bar " + std::to_string(b) + " is already active — activation must be unambiguous");
      for (Id n : add_nodes[pi])
        if (!active_sup.insert(n).second)
          errs.push_back(P + "support node " + std::to_string(n) + " is already active — activation must be unambiguous");
      for (Id b : ph.deactivate_bars) {
        if (!active_bars.erase(b)) {
          errs.push_back(P + "cannot deactivate bar " + std::to_string(b) + " — it is not active");
          continue;
        }
        for (const Tendon& t : model->tendons)
          if (live.count(t.id) && std::count(t.bars.begin(), t.bars.end(), b))
            errs.push_back(P + "cannot deactivate bar " + std::to_string(b) + " — it carries the live tendon " +
                           std::to_string(t.id));
      }
      for (Id n : ph.remove_support_nodes)
        if (!active_sup.erase(n))
          errs.push_back(P + "cannot deactivate support node " + std::to_string(n) + " — it is not active");
      for (Id c : ph.cases) {
        const auto lc = std::find_if(model->cases.begin(), model->cases.end(), [&](const LoadCase& x) { return x.id == c; });
        if (lc == model->cases.end()) {
          errs.push_back(P + "unknown load case " + std::to_string(c));
          continue;
        }
        if (tendon_case.count(c))
          errs.push_back(P + "case " + std::to_string(c) +
                         " contains a tendon load — tendon loads enter via tension orders, not cases");
        for (const Load& l : lc->loads)
          if (l.bar != 0 && !active_bars.count(l.bar)) {
            errs.push_back(P + "case " + std::to_string(c) + " applies a load to inactive bar " + std::to_string(l.bar));
            break;
          }
      }
      for (const TensionAction& ta : ph.tension) {
        const Tendon& t = model->tendon(ta.tendon);
        for (Id b : t.bars)
          if (!active_bars.count(b)) {
            errs.push_back(P + "cannot tension tendon " + std::to_string(t.id) + " — its chain bar " +
                           std::to_string(b) + " is not active");
            break;
          }
        if (ta.fraction > 0) live.insert(t.id);
        else live.erase(t.id);
      }
    }
    if (!errs.empty()) {
      std::string msg = "phase plan invalid: ";
      for (std::size_t i = 0; i < errs.size(); ++i) msg += (i ? "; " : "") + errs[i];
      throw ComponentError(msg);
    }
  }

  if (ctx.cancelled()) throw ComponentError::cancelled();
  PhasingOptions opt;
  opt.stations_per_bar = stations;
  const SolveHooks hooks = ctx_hooks(ctx, info().name);
  opt.hooks = &hooks;
  auto out = std::make_shared<PhasedModel>();
  out->model = model;
  try {
    std::lock_guard<std::mutex> lk(cache->mu);
    out->phased = solve_phased_cached(work, plan, opt, cache->cache);
  } catch (const Cancelled&) {
    throw ComponentError::cancelled();
  } catch (const std::exception& e) {
    throw ComponentError(std::string("phased analysis failed: ") + e.what());
  }
  // The two views: cases = phases (id = 1-based index, name = phase name).
  auto view = [&](bool cumulated) {
    auto a = std::make_shared<AnalyzedModel>();
    a->model = model;
    for (std::size_t i = 0; i < out->phased.phases.size(); ++i) {
      const PhaseResult& pr = out->phased.phases[i];
      CaseResult c = cumulated ? pr.cumulated : pr.incremental;
      c.id = static_cast<Id>(i + 1);
      c.name = pr.name;
      a->results.cases.push_back(std::move(c));
    }
    // Section-point stresses of the phase view, on demand (one shared
    // section cache: the second view re-solves nothing).
    a->stress_source = std::make_shared<const LazySectionStress>(stress_cache);
    return ResultsRef(std::move(a));
  };
  out->cumulated = view(true);
  out->incremental = view(false);
  if (ctx.progress) ctx.progress->report(info().name, "done", 1.0);
  ResultsRef cum = out->cumulated, inc = out->incremental;
  return {single(Value::of(cum)), single(Value::of(inc)), single(Value::of(PhasedRef(std::move(out))))};
}
json::Value PhasedAnalyzeNode::params() const {
  json::Value p = json::Value::object();
  p.set("stations", stations);
  p.set("tolerance", tolerance);
  return p;
}
void PhasedAnalyzeNode::set_params(const json::Value& p) {
  if (const auto* s = p.get("stations"); s && s->is_number() && s->as_number() >= 2 &&
                                         s->as_number() == std::trunc(s->as_number()))
    stations = static_cast<int>(s->as_number());
  if (const auto* t = p.get("tolerance"); t && t->is_number()) tolerance = t->as_number();
}

}  // namespace keta::param
