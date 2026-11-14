#include "keta/param/value.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace keta::param {

namespace {

bool same_concrete(const std::optional<ConcreteTime>& a, const std::optional<ConcreteTime>& b) {
  if (a.has_value() != b.has_value()) return false;
  if (!a) return true;
  return a->fck == b->fck && a->rh == b->rh && a->h0 == b->h0 && a->cement_s == b->cement_s &&
         a->t_shrink0 == b->t_shrink0;
}

bool same_props(const SectionProps& a, const SectionProps& b) {
  return a.A == b.A && a.Avy == b.Avy && a.Avz == b.Avz && a.Iy == b.Iy && a.Iz == b.Iz &&
         a.J == b.J && a.Wely == b.Wely && a.Welz == b.Welz && a.Wply == b.Wply && a.Wplz == b.Wplz;
}

bool same_shape(const Shape& a, const Shape& b) {
  if (a.index() != b.index()) return false;
  return std::visit(
      [&](const auto& x) -> bool {
        using T = std::decay_t<decltype(x)>;
        const T& y = std::get<T>(b);
        if constexpr (std::is_same_v<T, ShapeGeneric>) {
          return true;
        } else if constexpr (std::is_same_v<T, ShapeRect>) {
          return x.b == y.b && x.h == y.h;
        } else if constexpr (std::is_same_v<T, ShapeCircle>) {
          return x.d == y.d;
        } else if constexpr (std::is_same_v<T, ShapeTube>) {
          return x.d == y.d && x.t == y.t;
        } else if constexpr (std::is_same_v<T, ShapePolygon>) {
          return x.outer == y.outer && x.holes == y.holes;
        } else if constexpr (std::is_same_v<T, ShapeI>) {
          return x.h == y.h && x.b == y.b && x.tw == y.tw && x.tf == y.tf && x.r == y.r;
        } else {  // ShapeBox, ShapeT share {h, b, tw, tf}
          return x.h == y.h && x.b == y.b && x.tw == y.tw && x.tf == y.tf;
        }
      },
      a);
}

bool same_points(const std::vector<StressPoint>& a, const std::vector<StressPoint>& b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i)
    if (a[i].name != b[i].name || a[i].y != b[i].y || a[i].z != b[i].z ||
        a[i].tau_factor_y != b[i].tau_factor_y || a[i].tau_factor_z != b[i].tau_factor_z)
      return false;
  return true;
}

std::string fmt(double x) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%g", x);
  return buf;
}
std::string fmt3(const Vec3& v) { return "(" + fmt(v.x()) + ", " + fmt(v.y()) + ", " + fmt(v.z()) + ")"; }

}  // namespace

SetSpec SetSpec::make(std::string name, std::vector<std::int64_t> idx) {
  std::sort(idx.begin(), idx.end());
  idx.erase(std::unique(idx.begin(), idx.end()), idx.end());
  return {std::move(name), std::move(idx)};
}

bool materials_equal(const Material& a, const Material& b) {
  return a.name == b.name && a.E == b.E && a.nu == b.nu && a.rho == b.rho && a.alpha == b.alpha &&
         a.fy == b.fy && a.damping == b.damping && same_concrete(a.concrete, b.concrete);
}

bool sections_equal(const Section& a, const Section& b) {
  return a.name == b.name && same_shape(a.shape, b.shape) && same_props(a.props, b.props) &&
         same_points(a.stress_points, b.stress_points) && a.exact == b.exact;
}

std::optional<double> Value::as_number() const {
  if (const auto* x = get<double>()) return *x;
  if (const auto* i = get<std::int64_t>()) return static_cast<double>(*i);
  return std::nullopt;
}

std::optional<std::int64_t> Value::as_int() const {
  if (const auto* i = get<std::int64_t>()) return *i;
  if (const auto* x = get<double>())
    if (std::isfinite(*x) && std::trunc(*x) == *x && std::fabs(*x) < 9.0e15)
      return static_cast<std::int64_t>(*x);
  return std::nullopt;
}

std::optional<bool> Value::as_bool() const {
  if (const auto* b = get<bool>()) return *b;
  if (const auto* x = get<double>()) return *x != 0.0;
  if (const auto* i = get<std::int64_t>()) return *i != 0;
  return std::nullopt;
}

std::optional<Vec3> Value::as_point() const {
  if (const auto* p = get<Point>()) return p->p;
  if (const auto* v = get<Vector>()) return v->v;
  return std::nullopt;
}

std::string Value::type_tag() const {
  struct V {
    std::string operator()(double) const { return "number"; }
    std::string operator()(std::int64_t) const { return "int"; }
    std::string operator()(bool) const { return "bool"; }
    std::string operator()(const std::string&) const { return "text"; }
    std::string operator()(const Point&) const { return "point"; }
    std::string operator()(const Vector&) const { return "vector"; }
    std::string operator()(const Line&) const { return "line"; }
    std::string operator()(const MaterialSpec&) const { return "material"; }
    std::string operator()(const SectionSpec&) const { return "section"; }
    std::string operator()(const BeamSpec&) const { return "beam"; }
    std::string operator()(const SupportSpec&) const { return "support"; }
    std::string operator()(const PointMassSpec&) const { return "pointmass"; }
    std::string operator()(const LoadSpec&) const { return "load"; }
    std::string operator()(const ComboSpec&) const { return "combo"; }
    std::string operator()(const ModelRef&) const { return "model"; }
    std::string operator()(const ResultsRef&) const { return "results"; }
    std::string operator()(const FoundationSpec&) const { return "foundation"; }
    std::string operator()(const SetSpec&) const { return "set"; }
    std::string operator()(const BucklingRef&) const { return "buckling"; }
    std::string operator()(const ModalRef&) const { return "modal"; }
    std::string operator()(const EnvelopeRef&) const { return "envelope"; }
    std::string operator()(const TendonSpec&) const { return "tendon"; }
    std::string operator()(const PhaseSpec&) const { return "phase"; }
    std::string operator()(const PhasedRef&) const { return "phased"; }
    std::string operator()(const PathSpec&) const { return "path"; }
    std::string operator()(const ConvoySpec&) const { return "convoy"; }
    std::string operator()(const InfluenceRef&) const { return "influence"; }
    std::string operator()(const CableSpec&) const { return "cable"; }
    std::string operator()(const CableStaticRef&) const { return "cablestatic"; }
    std::string operator()(const HarmonicRef&) const { return "harmonic"; }
    std::string operator()(const ArcTraceRef&) const { return "arctrace"; }
    std::string operator()(const SlabSpec&) const { return "slab"; }
    std::string operator()(const NonlinearRef&) const { return "nonlinear"; }
    std::string operator()(const PsdRef&) const { return "psd"; }
    std::string operator()(const SpectrumRef&) const { return "spectrum"; }
  };
  return std::visit(V{}, d_);
}

std::string Value::describe() const {
  struct V {
    std::string operator()(double x) const { return fmt(x); }
    std::string operator()(std::int64_t i) const { return std::to_string(i); }
    std::string operator()(bool b) const { return b ? "True" : "False"; }
    std::string operator()(const std::string& s) const { return s; }
    std::string operator()(const Point& p) const { return "Point " + fmt3(p.p); }
    std::string operator()(const Vector& v) const { return "Vector " + fmt3(v.v); }
    std::string operator()(const Line& l) const { return "Line " + fmt3(l.a) + " -> " + fmt3(l.b); }
    std::string operator()(const MaterialSpec& m) const {
      return "Material " + (m.m.name.empty() ? std::string("(unnamed)") : m.m.name) + " E=" + fmt(m.m.E);
    }
    std::string operator()(const SectionSpec& s) const { return "Section " + s.s.name + " A=" + fmt(s.s.props.A); }
    std::string operator()(const BeamSpec& b) const {
      return std::string(b.truss ? "Truss " : "Beam ") + fmt3(b.line.a) + " -> " + fmt3(b.line.b);
    }
    std::string operator()(const SupportSpec& s) const {
      return s.line ? "Line support " + fmt3(s.line->a) + " -> " + fmt3(s.line->b)
                    : "Support " + fmt3(s.point);
    }
    std::string operator()(const PointMassSpec& m) const { return "Point mass " + fmt(m.mass) + " kg"; }
    std::string operator()(const LoadSpec& l) const { return "Load (case " + std::to_string(l.case_id) + ")"; }
    std::string operator()(const ComboSpec& c) const {
      return "Combination " + c.name.value_or("") + " (" + std::to_string(c.terms.size()) + " terms)";
    }
    std::string operator()(const ModelRef& m) const {
      return m ? "Model: " + std::to_string(m->nodes.size()) + " nodes, " + std::to_string(m->bars.size()) +
                     " bars, " + std::to_string(m->cases.size()) + " cases"
               : "Model (null)";
    }
    std::string operator()(const ResultsRef& r) const {
      return r ? "Results: " + std::to_string(r->results.cases.size()) + " cases, " +
                     std::to_string(r->results.combos.size()) + " combinations"
               : "Results (null)";
    }
    std::string operator()(const FoundationSpec& f) const {
      return "Foundation on " + std::to_string(f.lines.size()) + " beams";
    }
    std::string operator()(const SetSpec& s) const {
      return "Set " + s.name + " (" + std::to_string(s.indices.size()) + " members)";
    }
    std::string operator()(const BucklingRef& b) const {
      return b ? "Buckling: " + std::to_string(b->results.load_factors.size()) + " factors" : "Buckling (null)";
    }
    std::string operator()(const ModalRef& m) const {
      return m ? "Modal: " + std::to_string(m->results.frequencies_hz.size()) + " modes" : "Modal (null)";
    }
    std::string operator()(const TendonSpec& t) const {
      return "Tendon" + (t.name ? " " + *t.name : std::string()) + " on " + std::to_string(t.lines.size()) +
             (t.lines.size() == 1 ? " bar" : " bars") + ", case " + std::to_string(t.case_id);
    }
    std::string operator()(const PhaseSpec& p) const {
      return "Phase" + (p.name ? " " + *p.name : std::string()) + " (" + std::to_string(p.cases.size()) + " cases, " +
             std::to_string(p.tension.size()) + " tendons)";
    }
    std::string operator()(const PathSpec& p) const {
      return "Path" + (p.name.empty() ? std::string() : " " + p.name) + " (" + std::to_string(p.lines.size()) +
             (p.lines.size() == 1 ? " beam)" : " beams)");
    }
    std::string operator()(const ConvoySpec& c) const {
      return "Convoy" + (c.name.empty() ? std::string() : " " + c.name) + " (" + std::to_string(c.axles.size()) +
             " axles, " + std::to_string(c.segments.size()) + " segments)";
    }
    std::string operator()(const InfluenceRef& i) const {
      return i ? "Influence line " + i->label + " (" + std::to_string(i->line.s.size()) + " stations)"
               : "Influence (null)";
    }
    std::string operator()(const CableSpec& c) const { return "Cable " + fmt3(c.a) + " -> " + fmt3(c.b); }
    std::string operator()(const CableStaticRef& c) const {
      if (!c) return "Cable static (null)";
      std::size_t n = 0;
      for (const Bar& b : c->model->bars) n += b.kind == BarKind::Cable;
      return "Cable static: " + std::to_string(n) + (n == 1 ? " cable, " : " cables, ") +
             std::to_string(c->results.iterations) + " iterations";
    }
    std::string operator()(const HarmonicRef& h) const {
      return h ? "Harmonic: " + std::to_string(h->results.freqs_hz.size()) + " frequencies" : "Harmonic (null)";
    }
    std::string operator()(const SlabSpec& sl) const {
      return "Slab t=" + fmt(sl.thickness) + " (" + std::to_string(sl.outline.size()) + "-point outline, " +
             std::to_string(sl.holes.size()) + (sl.holes.size() == 1 ? " hole)" : " holes)");
    }
    std::string operator()(const NonlinearRef& n) const {
      return n ? "Nonlinear: " + std::to_string(n->results.iterations) + " Newton iterations" : "Nonlinear (null)";
    }
    std::string operator()(const PsdRef& p) const {
      return p ? "PSD: " + std::to_string(p->results.freqs_hz.size()) + " frequencies" : "PSD (null)";
    }
    std::string operator()(const SpectrumRef& sp) const {
      return sp ? "Response spectrum: " + std::to_string(sp->results.frequencies_hz.size()) + " modes, base shear " +
                      fmt(sp->results.base_shear) + " N"
                : "Spectrum (null)";
    }
    std::string operator()(const ArcTraceRef& t) const {
      return t ? "Arc-length trace: " + std::to_string(t->trace.points.size()) + " points" : "Arc trace (null)";
    }
    std::string operator()(const PhasedRef& p) const {
      return p ? "Phased analysis: " + std::to_string(p->phased.phases.size()) + " phases" : "Phased (null)";
    }
    std::string operator()(const EnvelopeRef& e) const {
      return e ? "Envelope" + (e->name ? " " + *e->name : std::string()) + " over " +
                     std::to_string(e->sources.size()) + (e->sources.size() == 1 ? " source" : " sources")
               : "Envelope (null)";
    }
  };
  return std::visit(V{}, d_);
}

}  // namespace keta::param
