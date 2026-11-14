#include "keta/model.h"

#include "keta/section_exact.h"

#include <cmath>
#include <set>
#include <unordered_set>

namespace keta {

// ---------------- Load constructors ----------------

Load Load::nodal(Id node, const Vec3& F, const Vec3& M) {
  Load l;
  l.kind = Kind::NodalForce;
  l.node = node;
  l.F = F;
  l.M = M;
  return l;
}
Load Load::udl(Id bar, const Vec3& q, LoadFrame fr, double a, double b) {
  Load l;
  l.kind = Kind::BarUdl;
  l.bar = bar;
  l.F = q;
  l.frame = fr;
  l.a = a;
  l.b = b;
  return l;
}
Load Load::parabolic(Id bar, const Vec3& qa, const Vec3& qmid, const Vec3& qb, LoadFrame fr, double a,
                     double b) {
  Load l = Load::linear(bar, qa, qb, fr, a, b);
  l.kind = Kind::BarParabolic;
  l.Fmid = qmid;
  return l;
}

Load Load::linear(Id bar, const Vec3& qa, const Vec3& qb, LoadFrame fr, double a, double b) {
  Load l;
  l.kind = Kind::BarLinear;
  l.bar = bar;
  l.F = qa;
  l.F2 = qb;
  l.frame = fr;
  l.a = a;
  l.b = b;
  return l;
}
Load Load::point(Id bar, double a, const Vec3& F, const Vec3& M, LoadFrame fr) {
  Load l;
  l.kind = Kind::BarPoint;
  l.bar = bar;
  l.a = a;
  l.F = F;
  l.M = M;
  l.frame = fr;
  return l;
}
Load Load::thermal(Id bar, double dt, double gy, double gz) {
  Load l;
  l.kind = Kind::BarThermal;
  l.bar = bar;
  l.dt = dt;
  l.gy = gy;
  l.gz = gz;
  return l;
}
Load Load::self_weight(const Vec3& g) {
  Load l;
  l.kind = Kind::SelfWeight;
  l.g = g;
  return l;
}
Load Load::support_disp(Id node, const std::array<std::optional<double>, 6>& d) {
  Load l;
  l.kind = Kind::SupportDisp;
  l.node = node;
  l.disp = d;
  return l;
}
Load Load::surface_pressure(Id slab, const Vec3& w, SurfaceFrame fr) {
  Load l;
  l.kind = Kind::SurfacePressure;
  l.slab = slab;
  l.F = w;
  l.sframe = fr;
  return l;
}
Load Load::surface_thermal(Id slab, double dt, double grad) {
  Load l;
  l.kind = Kind::SurfaceThermal;
  l.slab = slab;
  l.dt = dt;
  l.grad = grad;
  return l;
}
Load Load::prestrain(Id bar, double strain, double kappa_y, double kappa_z) {
  Load l;
  l.kind = Kind::Prestrain;
  l.bar = bar;
  l.strain = strain;
  l.kappa_y = kappa_y;
  l.kappa_z = kappa_z;
  return l;
}
Load Load::prestrain_point(Id bar, double a, double dux, double theta_y, double theta_z) {
  Load l;
  l.kind = Kind::PrestrainPoint;
  l.frame = LoadFrame::Local;  // Global: the jumps are given in global axes
  l.bar = bar;
  l.a = a;
  l.strain = dux;
  l.kappa_y = theta_y;
  l.kappa_z = theta_z;
  return l;
}
Load Load::prestrain_linear(Id bar, double strain_i, double kappa_y_i, double kappa_z_i,
                            double strain_j, double kappa_y_j, double kappa_z_j) {
  Load l;
  l.kind = Kind::PrestrainLinear;
  l.bar = bar;
  l.strain = strain_i;
  l.kappa_y = kappa_y_i;
  l.kappa_z = kappa_z_i;
  l.strain2 = strain_j;
  l.kappa_y2 = kappa_y_j;
  l.kappa_z2 = kappa_z_j;
  return l;
}
Load Load::prestrain_parabolic(Id bar, const Vec3& at_i, const Vec3& at_mid, const Vec3& at_j) {
  Load l = prestrain_linear(bar, at_i(0), at_i(1), at_i(2), at_j(0), at_j(1), at_j(2));
  l.kind = Kind::PrestrainParabolic;
  l.strain_m = at_mid(0);
  l.kappa_y_m = at_mid(1);
  l.kappa_z_m = at_mid(2);
  return l;
}

// ---------------- builders ----------------

Id Model::add_node(double x, double y, double z) {
  Id id = nodes.empty() ? 1 : nodes.back().id + 1;
  nodes.push_back({id, Vec3(x, y, z)});
  indexed_ = false;
  return id;
}
Id Model::add_material(const Material& m) {
  Material copy = m;
  if (copy.id == 0) copy.id = materials.empty() ? 1 : materials.back().id + 1;
  materials.push_back(copy);
  indexed_ = false;
  return copy.id;
}
Id Model::add_section(const Section& s) {
  Section copy = s;
  if (copy.id == 0) copy.id = sections.empty() ? 1 : sections.back().id + 1;
  copy.finalize();
  sections.push_back(copy);
  indexed_ = false;
  return copy.id;
}
Id Model::add_bar(Id ni, Id nj, Id mat, Id sec, BarKind kind) {
  Bar b;
  b.id = bars.empty() ? 1 : bars.back().id + 1;
  b.ni = ni;
  b.nj = nj;
  b.material = mat;
  b.section = sec;
  b.kind = kind;
  bars.push_back(b);
  indexed_ = false;
  resolve_shear(bars.back());
  return b.id;
}

void Model::resolve_shear(Bar& b) {
  if (b.kind != BarKind::Frame) return;
  Section* s = nullptr;
  for (Section& x : sections)
    if (x.id == b.section) s = &x;
  const Material* mt = nullptr;
  for (const Material& x : materials)
    if (x.id == b.material) mt = &x;
  if (!s || !mt || std::holds_alternative<ShapeGeneric>(s->shape) || s->exact) return;
  if (default_shear_areas() == ShearAreas::ShapeFactor) return;
  const double nu = mt->nu;
  if (std::isnan(s->shear_nu)) {
    // the closed-form placeholders still in place: the energy areas for this
    // material (a hand-set value is the user's)
    if (s->props.Avy != s->avy_shape || s->props.Avz != s->avz_shape) return;
    double ay = 0, az = 0;
    if (!energy_shear_areas(s->shape, nu, ay, az)) return;
    s->props.Avy = ay;
    s->props.Avz = az;
    s->shear_nu = nu;
    return;
  }
  if (s->shear_nu == nu) return;
  // another Poisson ratio on the same section: its own copy
  const Id base = s->shear_base ? s->shear_base : s->id;
  for (const Section& x : sections)
    if ((x.shear_base == base || x.id == base) && x.shear_nu == nu) {
      b.section = x.id;
      return;
    }
  Section c = *s;
  c.id = sections.back().id + 1;
  c.shear_base = base;
  double ay = 0, az = 0;
  energy_shear_areas(c.shape, nu, ay, az);
  c.props.Avy = ay;
  c.props.Avz = az;
  c.shear_nu = nu;
  sections.push_back(c);
  b.section = c.id;
}
void Model::add_support(Id node, const std::array<bool, 6>& fixed,
                        const std::array<double, 6>& spring) {
  Support s;
  s.node = node;
  s.fixed = fixed;
  s.spring = spring;
  supports.push_back(s);
}
Id Model::add_case(const std::string& name, std::vector<Load> loads) {
  Id id = cases.empty() ? 1 : cases.back().id + 1;
  cases.push_back({id, name, std::move(loads)});
  indexed_ = false;
  return id;
}
Id Model::add_combo(const std::string& name, std::vector<std::pair<Id, double>> terms) {
  Id id = combos.empty() ? 1 : combos.back().id + 1;
  combos.push_back({id, name, std::move(terms)});
  indexed_ = false;
  return id;
}
Id Model::add_slab(const Slab& s) {
  Slab copy = s;
  if (copy.id == 0) copy.id = slabs.empty() ? 1 : slabs.back().id + 1;
  slabs.push_back(std::move(copy));
  indexed_ = false;
  return slabs.back().id;
}
Id Model::add_shell(const std::array<Id, 3>& n, double t, Id mat, Id slab) {
  ShellTri sh;
  sh.id = shells.empty() ? 1 : shells.back().id + 1;
  sh.n = n;
  sh.t = t;
  sh.material = mat;
  sh.slab = slab;
  shells.push_back(sh);
  indexed_ = false;
  return sh.id;
}

// ---------------- link laws ----------------
namespace {
std::size_t law_segment(const std::vector<std::pair<double, double>>& p, double d) {
  // the segment [i, i+1] holding d (the end ones extended)
  std::size_t i = 0;
  while (i + 2 < p.size() && d >= p[i + 1].first) ++i;
  return i;
}
}  // namespace
double LinkLaw::tangent(double d) const {
  const std::size_t i = law_segment(pts, d);
  return (pts[i + 1].second - pts[i].second) / (pts[i + 1].first - pts[i].first);
}
double LinkLaw::force(double d) const {
  const std::size_t i = law_segment(pts, d);
  return pts[i].second + tangent(d) * (d - pts[i].first);
}
double LinkLaw::energy(double d) const {
  // int_0^d F: exact over the linear pieces between 0 and d
  auto prim = [&](double a, double b) {  // int_a^b F, a <= b
    double s = 0;
    std::vector<double> cuts{a};
    for (const auto& q : pts)
      if (q.first > a && q.first < b) cuts.push_back(q.first);
    cuts.push_back(b);
    for (std::size_t k = 0; k + 1 < cuts.size(); ++k)
      s += 0.5 * (force(cuts[k]) + force(cuts[k + 1])) * (cuts[k + 1] - cuts[k]);
    return s;
  };
  return d >= 0 ? prim(0.0, d) : -prim(d, 0.0);
}
LinkLaw LinkLaw::tension_only(double k) { return LinkLaw{{{-1.0, 0.0}, {0.0, 0.0}, {1.0, k}}}; }
LinkLaw LinkLaw::compression_only(double k) { return LinkLaw{{{-1.0, -k}, {0.0, 0.0}, {1.0, 0.0}}}; }
LinkLaw LinkLaw::gap(double g, double k) {
  return LinkLaw{{{-g - 1.0, -k}, {-g, 0.0}, {g, 0.0}, {g + 1.0, k}}};
}
LinkLaw LinkLaw::bilinear(double k, double Fy, double kp) {
  const double dy = Fy / k;
  return LinkLaw{{{-dy - 1.0, -Fy - kp}, {-dy, -Fy}, {dy, Fy}, {dy + 1.0, Fy + kp}}};
}

Load Load::shell_prestress(Id tendon, double factor) {
  Load l;
  l.kind = Kind::ShellPrestress;
  l.node = tendon;
  l.strain = factor;
  return l;
}

Load Load::shell_imposed(Id shell, const Eigen::Matrix<double, 18, 1>& d) {
  Load l;
  l.kind = Kind::ShellImposed;
  l.node = shell;
  for (int k = 0; k < 6; ++k) l.poly.push_back(d.segment<3>(3 * k));
  return l;
}

Id Model::add_continuum(Continuum::Kind kind, const std::vector<Id>& nodes, Id material, double t) {
  Continuum c;
  Id mx = 0;
  for (const Continuum& x : continua) mx = std::max(mx, x.id);
  c.id = mx + 1;
  c.kind = kind;
  c.n = nodes;
  c.material = material;
  c.t = t;
  const std::size_t n = nodes.size();
  const bool solid = kind == Continuum::Kind::Solid;
  require(solid ? (n == 4 || n == 8 || n == 20) : (n == 3 || n == 4 || n == 6 || n == 8),
          "continuum: a " + std::string(solid ? "solid takes 4, 8 or 20" : "plane element takes 3, 4, 6 or 8") + " nodes");
  require(std::isfinite(t) && t > 0, "continuum: thickness must be > 0");
  (void)this->material(material);
  for (Id id : nodes) (void)node(id);
  continua.push_back(c);
  return c.id;
}

Load Load::hydrostatic(double gamma, const Vec3& surface_point, Id slab, double side, const Vec3& down) {
  Load l;
  l.kind = Kind::Hydrostatic;
  l.hydro_gamma = gamma;
  l.hydro_point = surface_point;
  l.slab = slab;
  l.hydro_side = side;
  l.g = down;
  return l;
}

Load Load::polygon_pressure(const std::vector<Vec3>& poly, const Vec3& q,
                            const std::vector<std::vector<Vec3>>& holes, Id slab, const Vec3& proj) {
  Load l;
  l.kind = Kind::PolygonPressure;
  l.poly = poly;
  l.holes = holes;
  l.F = q;
  l.slab = slab;
  l.proj = proj;
  return l;
}

Id Model::add_rigid_link(Id master, std::vector<Id> slaves, const std::array<bool, 6>& dofs) {
  RigidLink l;
  l.id = rigid_links.empty() ? 1 : rigid_links.back().id + 1;
  l.master = master;
  l.slaves = std::move(slaves);
  l.dofs = dofs;
  rigid_links.push_back(std::move(l));
  return rigid_links.back().id;
}
Id Model::add_constraint(std::vector<ConstraintTerm> terms) {
  Constraint c;
  c.id = constraints.empty() ? 1 : constraints.back().id + 1;
  c.terms = std::move(terms);
  constraints.push_back(std::move(c));
  invalidate_index();
  return constraints.back().id;
}
Id Model::add_skew_support(Id node, const Mat3& axes, const std::array<bool, 6>& fixed) {
  SkewSupport s;
  s.id = skew_supports.empty() ? 1 : skew_supports.back().id + 1;
  s.node = node;
  s.axes = axes;
  s.fixed = fixed;
  skew_supports.push_back(s);
  invalidate_index();
  return skew_supports.back().id;
}
Id Model::add_elastic_link(Id ni, Id nj, const Mat6& k, double s) {
  ElasticLink l;
  l.id = links.empty() ? 1 : links.back().id + 1;
  l.ni = ni;
  l.nj = nj;
  l.k = k;
  l.s = s;
  links.push_back(l);
  return l.id;
}

// ---------------- lookups ----------------

void Model::index() const {
  // The memo is trusted only while every map still spans its whole vector.
  // Drivers fork a model and edit the PUBLIC entity vectors in place — phasing
  // rebuilds `bars` per phase, firechain appends the hot sections/materials,
  // the ST1 EXECs replace `cases` / `tendons` — and a plain copy inherits the
  // parent's `indexed_` flag together with the parent's maps. Without this
  // size check those lookups answer from the STALE layout ("unknown load case
  // id" on an id that is right there in the vector). Seven integer compares per
  // lookup; the O(n) rebuild still happens only when something actually moved.
  if (indexed_ && node_ix_.size() == nodes.size() && mat_ix_.size() == materials.size() &&
      sec_ix_.size() == sections.size() && bar_ix_.size() == bars.size() &&
      case_ix_.size() == cases.size() && combo_ix_.size() == combos.size() &&
      tendon_ix_.size() == tendons.size() && slab_ix_.size() == slabs.size() &&
      shell_ix_.size() == shells.size())
    return;
  auto build = [](auto& map, const auto& vec, const char* what) {
    map.clear();
    map.reserve(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i) {
      auto [it, inserted] = map.emplace(vec[i].id, static_cast<int>(i));
      if (!inserted) fail(std::string("model: duplicate ") + what + " id " + std::to_string(vec[i].id));
    }
  };
  build(node_ix_, nodes, "node");
  build(mat_ix_, materials, "material");
  build(sec_ix_, sections, "section");
  build(bar_ix_, bars, "bar");
  build(case_ix_, cases, "load case");
  build(combo_ix_, combos, "combination");
  build(tendon_ix_, tendons, "tendon");
  build(slab_ix_, slabs, "slab");
  build(shell_ix_, shells, "shell");
  indexed_ = true;
}

namespace {
template <class Map, class Vec>
const auto& lookup(const Map& map, const Vec& vec, Id id, const char* what) {
  auto it = map.find(id);
  if (it == map.end()) fail(std::string("model: unknown ") + what + " id " + std::to_string(id));
  return vec[static_cast<std::size_t>(it->second)];
}
}  // namespace

const Node& Model::node(Id id) const { index(); return lookup(node_ix_, nodes, id, "node"); }
const Material& Model::material(Id id) const { index(); return lookup(mat_ix_, materials, id, "material"); }
const Section& Model::section(Id id) const { index(); return lookup(sec_ix_, sections, id, "section"); }
const Bar& Model::bar(Id id) const { index(); return lookup(bar_ix_, bars, id, "bar"); }
const LoadCase& Model::load_case(Id id) const { index(); return lookup(case_ix_, cases, id, "load case"); }
const Combination& Model::combo(Id id) const { index(); return lookup(combo_ix_, combos, id, "combination"); }
const Tendon& Model::tendon(Id id) const { index(); return lookup(tendon_ix_, tendons, id, "tendon"); }
const Slab& Model::slab(Id id) const { index(); return lookup(slab_ix_, slabs, id, "slab"); }

bool Model::has_surface(Id id) const {
  index();
  if (slab_ix_.count(id)) return true;
  for (const ShellTri& sh : shells)
    if (sh.slab == id) return true;
  return false;
}

int Model::node_index(Id id) const {
  index();
  auto it = node_ix_.find(id);
  if (it == node_ix_.end()) fail("model: unknown node id " + std::to_string(id));
  return it->second;
}
int Model::bar_index(Id id) const {
  index();
  auto it = bar_ix_.find(id);
  if (it == bar_ix_.end()) fail("model: unknown bar id " + std::to_string(id));
  return it->second;
}

double Model::bar_chord_length(const Bar& b) const {
  return (node(b.nj).x - node(b.ni).x).norm();
}

// ---------------- validation ----------------

void Model::validate() const {
  index();
  require(!nodes.empty(), "model: no nodes");
  for (const Node& n : nodes)
    require(finite(n.x), "model: node " + std::to_string(n.id) + " has non-finite coordinates");

  for (const Material& m : materials) {
    require(finite(m.E) && m.E > 0, "material '" + m.name + "': E must be positive");
    require(finite(m.nu) && m.nu > -1 && m.nu < 0.5, "material '" + m.name + "': nu out of range");
    require(finite(m.rho) && m.rho >= 0, "material '" + m.name + "': rho must be >= 0");
    require(finite(m.alpha) && finite(m.fy), "material '" + m.name + "': non-finite parameter");
    require(finite(m.G_given) && m.G_given >= 0, "material '" + m.name + "': G must be >= 0 (0 = isotropic)");
  }
  // An independent shear modulus is a bar property: the shell and continuum
  // laws are isotropic in (E, nu) and would silently ignore it.
  auto isotropic_only = [&](Id material, const std::string& what) {
    const Material& m = this->material(material);
    require(m.G_given == 0, what + " uses material '" + m.name +
                                "' with an independent G: shells and continua take an isotropic (E, nu) material");
  };
  for (const ShellTri& sh : shells) isotropic_only(sh.material, "shell " + std::to_string(sh.id));
  for (const Continuum& c : continua) isotropic_only(c.material, "continuum element " + std::to_string(c.id));
  for (const Section& s : sections) {
    const SectionProps& p = s.props;
    require(finite(p.A) && p.A > 0, "section '" + s.name + "': A must be positive");
    for (double v : {p.Iy, p.Iz, p.J, p.Avy, p.Avz})
      require(finite(v) && v >= 0, "section '" + s.name + "': negative or non-finite property");
  }

  for (const ShellTendon& t : shell_tendons) {
    const std::string w = "shell tendon " + std::to_string(t.id) + ": ";
    require(t.pts.size() >= 2, w + "a tendon needs two points at least");
    for (std::size_t i = 0; i < t.pts.size(); ++i) {
      require(finite(t.pts[i]), w + "non-finite point");
      if (i > 0) require((t.pts[i] - t.pts[i - 1]).norm() > 0, w + "repeated point");
    }
    require(finite(t.P0) && t.P0 >= 0 && finite(t.mu) && t.mu >= 0 && finite(t.k) && t.k >= 0,
            w + "force and friction must be finite and >= 0");
    require(t.slab == 0 || has_surface(t.slab), w + "its slab matches no shell");
  }
  for (const Constraint& c : constraints) {
    const std::string w = "constraint " + std::to_string(c.id) + ": ";
    require(!c.terms.empty(), w + "no term");
    std::set<std::pair<Id, int>> seen;
    for (const ConstraintTerm& t : c.terms) {
      (void)node(t.node);
      require(t.comp >= 0 && t.comp < 6, w + "component 0..5");
      require(finite(t.coef), w + "non-finite coefficient");
      require(seen.insert({t.node, t.comp}).second, w + "a DOF listed twice");
    }
    require(c.terms.front().coef != 0, w + "the first (eliminated) term needs a non-zero coefficient");
  }
  for (const SkewSupport& sk : skew_supports) {
    const std::string w = "skew support " + std::to_string(sk.id) + ": ";
    (void)node(sk.node);
    require(sk.axes.allFinite() && (sk.axes * sk.axes.transpose() - Mat3::Identity()).cwiseAbs().maxCoeff() < 1e-9,
            w + "axes must be orthonormal (rows = local x, y, z)");
    bool any = false;
    for (bool b : sk.fixed) any = any || b;
    require(any, w + "no blocked component");
  }
  for (const RigidLink& l : rigid_links) {
    const std::string w = "rigid link " + std::to_string(l.id) + ": ";
    (void)node(l.master);
    require(!l.slaves.empty(), w + "no slave node");
    bool any = false;
    for (bool b : l.dofs) any = any || b;
    require(any, w + "no coupled component");
    std::set<Id> seen;
    for (Id s : l.slaves) {
      (void)node(s);
      require(s != l.master, w + "node " + std::to_string(s) + " is its own master");
      require(seen.insert(s).second, w + "slave " + std::to_string(s) + " listed twice");
    }
  }
  for (const ElasticLink& l : links) {
    const std::string w = "elastic link " + std::to_string(l.id) + ": ";
    (void)node(l.ni);
    (void)node(l.nj);
    require(l.ni != l.nj, w + "same node at both ends");
    require(l.k.allFinite() && (l.k - l.k.transpose()).cwiseAbs().maxCoeff() <= 1e-12 * l.k.cwiseAbs().maxCoeff(),
            w + "stiffness must be finite and symmetric");
    require(finite(l.s) && l.s >= 0 && l.s <= 1, w + "spring position s must lie in [0, 1]");
    require(finite(l.roll), w + "non-finite roll");
    for (int c = 0; c < 6; ++c) {
      const auto& law = l.law[static_cast<std::size_t>(c)];
      if (!law) continue;
      const auto& p = law->pts;
      require(p.size() >= 2, w + "a law needs two points at least");
      for (std::size_t i = 0; i < p.size(); ++i) {
        require(finite(p[i].first) && finite(p[i].second), w + "non-finite law point");
        if (i > 0) {
          require(p[i].first > p[i - 1].first, w + "law: deformations must increase strictly");
        }
      }
      for (int q = 0; q < 6; ++q)
        require(q == c || (l.k(c, q) == 0 && l.k(q, c) == 0), w + "a component with a law takes no coupling term in k");
    }
    if (l.axes) {
      const Mat3& R = *l.axes;
      require(R.allFinite() && (R * R.transpose() - Mat3::Identity()).cwiseAbs().maxCoeff() < 1e-9 &&
                  R.determinant() > 0,
              w + "axes must be a proper orthonormal frame (rows = local x, y, z)");
    }
  }

  for (const auto& [bid, cb] : composite) {
    const std::string w = "composite bar " + std::to_string(bid) + ": ";
    const Bar& b = bar(bid);
    require(b.kind == BarKind::Frame, w + "must be a frame bar");
    require(!b.has_ecc(), w + "no rigid arms (the parts carry their own positions)");
    require(!b.soil.any() && !tapered.count(bid), w + "no Winkler bed or exact tapered kernel");
    require(!cb.parts.empty(), w + "no part");
    for (const SectionPart& q : cb.parts) {
      (void)material(q.material);
      require(finite(q.props.A) && q.props.A > 0, w + "part '" + q.name + "': A must be positive");
      for (double v : {q.props.Iy, q.props.Iz, q.props.J, q.props.Avy, q.props.Avz, q.z})
        require(finite(v), w + "part '" + q.name + "': non-finite property");
      for (double v : {q.props.Iy, q.props.Iz, q.props.J, q.props.Avy, q.props.Avz})
        require(v >= 0, w + "part '" + q.name + "': negative property");
    }
  }

  std::set<std::pair<Id, Id>> seen_pairs;
  for (const Bar& b : bars) {
    const Node& ni = node(b.ni);
    const Node& nj = node(b.nj);
    require(b.ni != b.nj, "bar " + std::to_string(b.id) + ": zero-length (same node twice)");
    double L = (nj.x - ni.x).norm();
    require(L > 1e-12, "bar " + std::to_string(b.id) + ": degenerate length");
    (void)material(b.material);
    const Section& sec = section(b.section);
    if (b.kind == BarKind::Frame)
      require(sec.props.Iy > 0 && sec.props.Iz > 0,
              "bar " + std::to_string(b.id) + ": frame bar needs Iy, Iz > 0");
    require(finite(b.roll), "bar " + std::to_string(b.id) + ": non-finite roll");
    require(finite(b.ecc[0]) && finite(b.ecc[1]),
            "bar " + std::to_string(b.id) + ": non-finite eccentricity");
    for (double k : {b.soil.kfx, b.soil.kfy, b.soil.kfz, b.soil.kmx})
      require(finite(k) && k >= 0, "bar " + std::to_string(b.id) + ": negative/non-finite soil modulus");
    if (b.kind == BarKind::Truss) {
      require(!b.soil.any(), "bar " + std::to_string(b.id) + ": Winkler bed on a truss bar");
      require(!b.has_ecc(), "bar " + std::to_string(b.id) + ": eccentricity on a truss bar");
    }
    if (b.kind == BarKind::Cable) {
      require(b.cable.has_value(), "bar " + std::to_string(b.id) + ": cable bar without cable data");
      const CableProps& cp = *b.cable;
      require(std::isfinite(cp.ea) && cp.ea >= 0 && std::isfinite(cp.end_tension) && cp.end_tension >= 0 &&
                  (!cp.w || (std::isfinite(*cp.w) && *cp.w >= 0)),
              "bar " + std::to_string(b.id) + ": cable EA, w and end tension must be finite and non-negative");
    }
  }

  // ST1 CARA VAR: exact tapered member kernels (see `TaperedBar` in model.h).
  for (const auto& [bid, tb] : tapered) {
    const Bar& b = bar(bid);  // typed error on a dangling reference
    const std::string at = "tapered bar " + std::to_string(bid);
    require(b.kind == BarKind::Frame, at + ": a variable section needs a frame bar");
    for (bool r : b.releases)
      require(!r, at + ": end releases on a tapered bar are not supported (the release "
                       "condensation operator is built from the prismatic member block)");
    require(tb.ke.allFinite(), at + ": non-finite tapered stiffness");
  }

  for (const Support& s : supports) {
    (void)node(s.node);
    const std::string at = "support at node " + std::to_string(s.node);
    require(s.any(), at + ": no DOF fixed and no spring");
    for (double k : s.spring) require(finite(k) && k >= 0, at + ": bad spring");
    // ST1 APPUI EL full 6x6 + ROT alpha (see springmass.h).
    bool has_oneside = false;
    for (int c : s.oneside) has_oneside |= (c != 0);
    if (s.spring6) {
      require(s.spring6->allFinite(), at + ": non-finite full support stiffness");
      for (int c = 0; c < 6; ++c)
        require((*s.spring6)(c, c) >= 0, at + ": negative full-support diagonal term");
      require(!has_oneside,
              at + ": per-DOF one-sided release cannot be combined with a full 6x6 elastic "
                   "support (use DECOL, which releases the whole support)");
    }
    require(finite(s.rot), at + ": non-finite support rotation alpha");
    for (int c = 0; c < 6; ++c) {
      const double cp = s.cap[static_cast<std::size_t>(c)];
      require(finite(cp) && cp >= 0, at + ": support force cap must be finite and >= 0");
      require(cp == 0 || s.fixed[static_cast<std::size_t>(c)], at + ": a force cap needs a blocked component");
    }
    // ST1 DECOL: one global axis, all-or-nothing.
    require(s.decol_axis >= -1 && s.decol_axis <= 2, at + ": DECOL axis out of range (-1..2)");
    if (s.decol_axis >= 0) {
      require(s.decol_sign == 1 || s.decol_sign == -1, at + ": DECOL sign must be +1 (POS) or -1 (NEG)");
      require(!has_oneside, at + ": DECOL cannot be combined with per-DOF one-sided releases");
    }
  }

  for (const NodalMass& nm : nodal_masses) {
    (void)node(nm.node);
    const std::string at = "nodal mass at node " + std::to_string(nm.node);
    require(nm.m.allFinite(), at + ": non-finite mass matrix");
    for (int c = 0; c < 6; ++c) require(nm.m(c, c) >= 0, at + ": negative diagonal mass term");
  }
  for (const BarMass& bm : bar_masses) {
    (void)bar(bm.bar);
    const std::string at = "bar mass on bar " + std::to_string(bm.bar);
    for (double v : {bm.m1, bm.m2, bm.m3})
      require(finite(v) && v >= 0, at + ": negative/non-finite line mass");
    require(finite(bm.ecc_or) && finite(bm.ecc_ex), at + ": non-finite mass eccentricity");
  }

  // ---- flat shells: slab records ----
  for (const Slab& s : slabs) {
    const std::string at = "slab " + std::to_string(s.id);
    (void)material(s.material);
    require(finite(s.t) && s.t > 0, at + ": thickness must be finite and positive");
    require(finite(s.mesh_size) && s.mesh_size > 0, at + ": mesh size must be finite and positive");
    require(s.outline.size() >= 3, at + ": outline needs at least 3 points");
    for (const Vec3& p : s.outline) require(finite(p), at + ": non-finite outline point");
    for (const auto& h : s.holes) {
      require(h.size() >= 3, at + ": a hole ring needs at least 3 points");
      for (const Vec3& p : h) require(finite(p), at + ": non-finite hole point");
    }
    // Planarity: Newell normal of the outline, then the largest out-of-plane
    // deviation of any outline/hole point, relative to the bounding-box size.
    Vec3 nrm = Vec3::Zero();
    Vec3 lo = s.outline.front(), hi = s.outline.front();
    for (std::size_t i = 0; i < s.outline.size(); ++i) {
      const Vec3& a = s.outline[i];
      const Vec3& b = s.outline[(i + 1) % s.outline.size()];
      nrm += a.cross(b);
      lo = lo.cwiseMin(a);
      hi = hi.cwiseMax(a);
    }
    double diag = (hi - lo).norm();
    require(nrm.norm() > 1e-12 * std::max(diag * diag, 1e-12), at + ": degenerate (zero-area) outline");
    Vec3 unit_n = nrm.normalized();
    double dev = 0;
    auto sweep_plane = [&](const std::vector<Vec3>& ring) {
      for (const Vec3& p : ring) dev = std::max(dev, std::fabs((p - s.outline.front()).dot(unit_n)));
    };
    sweep_plane(s.outline);
    for (const auto& h : s.holes) sweep_plane(h);
    double tol = 1e-6 * std::max(diag, 1.0);
    require(dev <= tol, at + ": outline is not planar (deviation " + std::to_string(dev) +
                            " > tolerance " + std::to_string(tol) + ")");
    if (s.local_x) {
      require(finite(*s.local_x), at + ": non-finite local axis (local_x) override");
      double len = s.local_x->norm();
      require(len > 0, at + ": zero-length local axis (local_x) override");
      Vec3 au = *s.local_x / len;
      double in_plane = (au - unit_n * au.dot(unit_n)).norm();
      require(in_plane >= 1e-6, at + ": local axis (local_x) override is parallel to the slab normal");
    }
  }

  // ---- flat shells: lowered facets ----
  for (const ShellTri& sh : shells) {
    const std::string at = "shell " + std::to_string(sh.id);
    for (Id nid : sh.n) (void)node(nid);
    require(sh.n[0] != sh.n[1] && sh.n[1] != sh.n[2] && sh.n[0] != sh.n[2],
            at + ": repeated corner node");
    (void)material(sh.material);
    require(finite(sh.t) && sh.t > 0, at + ": thickness must be finite and positive");
    // The slab reference is enforced only in a slab-carrying model (a
    // standalone hand-built mesh uses `slab` as a free grouping id).
    if (!slabs.empty()) (void)slab(sh.slab);
  }

  // shell-tendon ids, for the ShellPrestress loads (built on first use)
  std::unordered_set<Id> shell_tendon_ids;
  bool shell_tendon_ids_built = false;
  for (const LoadCase& lc : cases) {
    for (const Load& l : lc.loads) {
      switch (l.kind) {
        case Load::Kind::NodalForce:
          (void)node(l.node);
          require(finite(l.F) && finite(l.M), "case '" + lc.name + "': non-finite nodal load");
          break;
        case Load::Kind::BarUdl:
        case Load::Kind::BarLinear:
        case Load::Kind::BarParabolic:
        case Load::Kind::BarPoint: {
          if (l.kind == Load::Kind::BarParabolic)
            require(l.M.isZero() && l.M2.isZero() && finite(l.Fmid),
                    "case '" + lc.name + "': a parabolic bar load carries forces only (finite, no couple)");
          const Bar& b = bar(l.bar);
          require(b.kind != BarKind::Cable, "case '" + lc.name + "': bar load on a cable");
          require(finite(l.F) && finite(l.M) && finite(l.F2),
                  "case '" + lc.name + "': non-finite bar load");
          require(finite(l.a) && finite(l.b), "case '" + lc.name + "': non-finite abscissa");
          break;
        }
        case Load::Kind::BarThermal: {
          const Bar& b = bar(l.bar);
          const Material& m = material(b.material);
          require(m.alpha != 0.0, "case '" + lc.name + "': thermal load on bar " +
                                      std::to_string(b.id) + " whose material has no alpha");
          require(finite(l.dt) && finite(l.gy) && finite(l.gz),
                  "case '" + lc.name + "': non-finite thermal load");
          break;
        }
        case Load::Kind::SelfWeight:
          require(finite(l.g), "case '" + lc.name + "': non-finite gravity vector");
          break;
        case Load::Kind::SupportDisp: {
          (void)node(l.node);
          bool any = false;
          for (const auto& d : l.disp)
            if (d) {
              any = true;
              require(finite(*d), "case '" + lc.name + "': non-finite prescribed displacement");
            }
          require(any, "case '" + lc.name + "': SupportDisp with no prescribed DOF");
          break;
        }
        case Load::Kind::Prestrain:
          (void)bar(l.bar);
          require(finite(l.strain) && finite(l.kappa_y) && finite(l.kappa_z),
                  "case '" + lc.name + "': non-finite prestrain");
          break;
        case Load::Kind::PrestrainPoint: {
          const Bar& b = bar(l.bar);
          require(b.kind != BarKind::Cable,
                  "case '" + lc.name + "': prestrain dislocation on a cable");
          require(finite(l.strain) && finite(l.kappa_y) && finite(l.kappa_z) && finite(l.a),
                  "case '" + lc.name + "': non-finite prestrain dislocation");
          break;
        }
        case Load::Kind::PrestrainLinear:
          (void)bar(l.bar);
          require(finite(l.strain) && finite(l.kappa_y) && finite(l.kappa_z) &&
                      finite(l.strain2) && finite(l.kappa_y2) && finite(l.kappa_z2),
                  "case '" + lc.name + "': non-finite linear prestrain");
          break;
        case Load::Kind::PrestrainParabolic:
          (void)bar(l.bar);
          require(finite(l.strain) && finite(l.kappa_y) && finite(l.kappa_z) && finite(l.kappa_x) &&
                      finite(l.strain2) && finite(l.kappa_y2) && finite(l.kappa_z2) && finite(l.kappa_x2) &&
                      finite(l.strain_m) && finite(l.kappa_y_m) && finite(l.kappa_z_m) && finite(l.kappa_x_m),
                  "case '" + lc.name + "': non-finite parabolic prestrain");
          break;
        case Load::Kind::SurfacePressure: {
          const std::string at = "case '" + lc.name + "': surface pressure on slab " +
                                 std::to_string(l.slab);
          require(has_surface(l.slab), at + " which matches no slab record and no shell group");
          require(finite(l.F), at + ": non-finite pressure");
          // A LocalNormal pressure carries its scalar in `w.z`; a non-zero
          // in-plane component would be silently dropped — reject it instead.
          if (l.sframe == SurfaceFrame::LocalNormal)
            require(l.F.x() == 0 && l.F.y() == 0,
                    at + ": a LocalNormal pressure carries its scalar in w.z; the in-plane "
                         "components must be zero");
          break;
        }
        case Load::Kind::ShellPrestress: {
          if (!shell_tendon_ids_built) {
            for (const ShellTendon& t : shell_tendons) shell_tendon_ids.insert(t.id);
            shell_tendon_ids_built = true;
          }
          const bool found = shell_tendon_ids.count(l.node) > 0;
          require(found, "case '" + lc.name + "': no shell tendon " + std::to_string(l.node));
          require(finite(l.strain), "case '" + lc.name + "': non-finite tendon factor");
          break;
        }
        case Load::Kind::ShellImposed: {
          const bool found = shell_ix_.count(l.node) > 0;  // indexed above (index())
          require(found, "case '" + lc.name + "': no shell " + std::to_string(l.node));
          require(l.poly.size() == 6, "case '" + lc.name + "': an imposed shell field holds 18 values");
          for (const Vec3& v : l.poly) require(finite(v), "case '" + lc.name + "': non-finite imposed shell field");
          break;
        }
        case Load::Kind::Hydrostatic: {
          const std::string at = "case '" + lc.name + "': hydrostatic pressure";
          require(l.slab == 0 || has_surface(l.slab), at + " on a slab that matches no shell");
          require(finite(l.hydro_gamma) && l.hydro_gamma >= 0, at + ": the unit weight must be finite and >= 0");
          require(finite(l.hydro_point) && finite(l.g) && l.g.norm() > 0, at + ": non-finite surface point or direction");
          require(l.hydro_side == 1 || l.hydro_side == -1, at + ": side must be +1 or -1");
          break;
        }
        case Load::Kind::PolygonPressure: {
          const std::string at = "case '" + lc.name + "': polygon pressure";
          require(l.slab == 0 || has_surface(l.slab), at + " on a slab that matches no shell");
          require(finite(l.F) && finite(l.proj) && l.proj.norm() > 0, at + ": non-finite load or direction");
          require(l.poly.size() >= 3, at + ": a polygon needs 3 points at least");
          for (const Vec3& p : l.poly) require(finite(p), at + ": non-finite vertex");
          for (const auto& h : l.holes) {
            require(h.size() >= 3, at + ": an opening needs 3 points at least");
            for (const Vec3& p : h) require(finite(p), at + ": non-finite opening vertex");
          }
          break;
        }
        case Load::Kind::SurfaceThermal: {
          const std::string at = "case '" + lc.name + "': surface thermal on slab " +
                                 std::to_string(l.slab);
          require(has_surface(l.slab), at + " which matches no slab record and no shell group");
          require(finite(l.dt) && finite(l.grad), at + ": non-finite temperature / gradient");
          // Resolve alpha through the Slab record when there is one, else through
          // the shells carrying the id. alpha == 0 is a typed error, never a
          // silent zero (the same strictness as the bar thermal load).
          bool checked = false;
          for (const Slab& s : slabs)
            if (s.id == l.slab) {
              require(material(s.material).alpha != 0.0,
                      at + " whose material has no alpha");
              checked = true;
            }
          if (!checked)
            for (const ShellTri& sh : shells)
              if (sh.slab == l.slab)
                require(material(sh.material).alpha != 0.0,
                        at + " whose material has no alpha");
          break;
        }
      }
    }
  }

  for (const Combination& c : combos)
    for (const auto& [cid, f] : c.terms) {
      (void)load_case(cid);
      require(finite(f), "combination '" + c.name + "': non-finite factor");
    }

  for (const Tendon& t : tendons) {
    require(!t.bars.empty(), "tendon '" + t.name + "': empty bar chain");
    for (Id bid : t.bars) {
      const Bar& b = bar(bid);
      require(b.kind == BarKind::Frame, "tendon '" + t.name + "': tendon on non-frame bar " +
                                            std::to_string(bid));
    }
    require(finite(t.jack_force) && t.jack_force >= 0, "tendon '" + t.name + "': bad jack force");
    require(!t.profile.empty(), "tendon '" + t.name + "': empty profile");
    require(t.flips.empty() || t.flips.size() == t.bars.size(),
            "tendon '" + t.name + "': flips must parallel the bar chain");
    for (std::size_t i = 0; i < t.force_profile.size(); ++i) {
      const auto& [fs, fp] = t.force_profile[i];
      require(finite(fs) && finite(fp) && fp >= 0, "tendon '" + t.name + "': bad force profile point");
      require(i == 0 || fs > t.force_profile[i - 1].first,
              "tendon '" + t.name + "': force profile abscissae must increase strictly");
    }
  }
}

}  // namespace keta
