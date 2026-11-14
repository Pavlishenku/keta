#include "keta/sections.h"

#include <algorithm>
#include <atomic>
#include <cmath>

namespace keta {

namespace {

// Torsion constant of a thin rectangle strip b*t (b >> t): J = b t^3 / 3 with
// the classical correction for finite aspect ratio when needed.
double j_rect(double b, double h) {
  // Saint-Venant's exact series for a solid rectangle, sides a >= c
  // (Timoshenko & Goodier, ch. 10):
  //   J = a c^3 / 3 [1 - 192/pi^5 (c/a) sum_(n odd) tanh(n pi a / 2c) / n^5]
  // (Roark's 1/3 - 0.21 r (1 - r^4/12) is its fit, 0.1 % off)
  const double a = std::max(b, h), c = std::min(b, h);
  double sum = 0;
  for (int n = 1; n < 400; n += 2) {
    const double term = std::tanh(n * kPi * a / (2 * c)) / std::pow(n, 5);
    sum += term;
    if (term < 1e-18 * sum) break;
  }
  return a * c * c * c / 3.0 * (1.0 - 192.0 / std::pow(kPi, 5) * (c / a) * sum);
}

void check_pos(double v, const char* what) {
  if (!(v > 0) || !std::isfinite(v)) fail(std::string("section: non-positive ") + what);
}

}  // namespace

SectionProps props_rect(const ShapeRect& s) {
  check_pos(s.b, "b");
  check_pos(s.h, "h");
  SectionProps p;
  p.A = s.b * s.h;
  p.Iz = s.h * s.b * s.b * s.b / 12.0;  // bending about z -> fiber distance along y
  p.Iy = s.b * s.h * s.h * s.h / 12.0;
  p.Avy = p.A * 5.0 / 6.0;
  p.Avz = p.A * 5.0 / 6.0;
  p.J = j_rect(s.b, s.h);
  p.Wely = p.Iy / (s.h / 2);
  p.Welz = p.Iz / (s.b / 2);
  p.Wply = s.b * s.h * s.h / 4.0;
  p.Wplz = s.h * s.b * s.b / 4.0;
  return p;
}

SectionProps props_circle(const ShapeCircle& s) {
  check_pos(s.d, "d");
  SectionProps p;
  double r = s.d / 2;
  p.A = kPi * r * r;
  p.Iy = p.Iz = kPi * r * r * r * r / 4.0;
  p.J = 2.0 * p.Iy;
  p.Avy = p.Avz = p.A * 0.9;
  p.Wely = p.Welz = p.Iy / r;
  p.Wply = p.Wplz = s.d * s.d * s.d / 6.0;
  return p;
}

SectionProps props_tube(const ShapeTube& s) {
  check_pos(s.d, "d");
  check_pos(s.t, "t");
  if (2 * s.t >= s.d) fail("section: tube wall thicker than radius");
  SectionProps p;
  double ro = s.d / 2, ri = ro - s.t;
  p.A = kPi * (ro * ro - ri * ri);
  p.Iy = p.Iz = kPi * (std::pow(ro, 4) - std::pow(ri, 4)) / 4.0;
  p.J = 2.0 * p.Iy;
  p.Avy = p.Avz = p.A / 2.0;  // thin tube shear area ~ A/2
  p.Wely = p.Welz = p.Iy / ro;
  double dm = s.d - s.t;
  p.Wply = p.Wplz = dm * dm * s.t;  // thin-wall plastic modulus d_m^2 t
  return p;
}

SectionProps props_i(const ShapeI& s) {
  check_pos(s.h, "h");
  check_pos(s.b, "b");
  check_pos(s.tw, "tw");
  check_pos(s.tf, "tf");
  if (2 * s.tf >= s.h) fail("section: I flanges overlap");
  SectionProps p;
  double hw = s.h - 2 * s.tf;  // clear web height
  p.A = 2 * s.b * s.tf + hw * s.tw;
  // Strong axis = y (bending in the plane of the web, fibers at +-h/2 along z).
  p.Iy = s.b * std::pow(s.h, 3) / 12.0 - (s.b - s.tw) * std::pow(hw, 3) / 12.0;
  p.Iz = 2 * s.tf * std::pow(s.b, 3) / 12.0 + hw * std::pow(s.tw, 3) / 12.0;
  p.J = (2 * s.b * std::pow(s.tf, 3) + hw * std::pow(s.tw, 3)) / 3.0;
  p.Avz = hw * s.tw;               // web carries Vz
  p.Avy = 2 * s.b * s.tf * 5. / 6.; // flanges carry Vy
  p.Wely = p.Iy / (s.h / 2);
  p.Welz = p.Iz / (s.b / 2);
  p.Wply = s.b * s.tf * (s.h - s.tf) + s.tw * hw * hw / 4.0;
  p.Wplz = s.tf * s.b * s.b / 2.0 + hw * s.tw * s.tw / 4.0;
  // warping constant of the doubly symmetric I on its flange mid-lines
  // (flanges as thin rectangles, web of no warping rigidity): the tabulated
  // EN 10365 value (IPE 300: 125.9e3 cm^6)
  p.Iw = s.tf * std::pow(s.b, 3) * std::pow(s.h - s.tf, 2) / 24.0;
  if (s.r == 0) return p;

  // Root fillets: four circular spandrels in the web-flange corners. Spandrel
  // of radius r (square r x r minus the quarter disk): area (1 - pi/4) r^2,
  // centroid at d = (5/6 - pi/4)/(1 - pi/4) r from each corner face, second
  // moment about a corner face (1 - 5 pi/16) r^4.
  check_pos(s.r, "r");
  if (2 * (s.tf + s.r) >= s.h || s.tw + 2 * s.r >= s.b) fail("section: I root fillets do not fit");
  const double r = s.r;
  const double As = (1.0 - kPi / 4.0) * r * r;
  const double d = (5.0 / 6.0 - kPi / 4.0) / (1.0 - kPi / 4.0) * r;
  const double I_own = (1.0 - 5.0 * kPi / 16.0) * std::pow(r, 4) - As * d * d;
  const double zc = s.h / 2 - s.tf - d;  // spandrel centroid, strong-axis lever
  const double yc = s.tw / 2 + d;        // weak-axis lever
  p.A += 4 * As;
  p.Iy += 4 * (I_own + As * zc * zc);
  p.Iz += 4 * (I_own + As * yc * yc);
  p.Wely = p.Iy / (s.h / 2);
  p.Welz = p.Iz / (s.b / 2);
  p.Wply += 4 * As * zc;
  p.Wplz += 4 * As * yc;
  // St-Venant torsion with the fillet bulbs (El Darwish & Johnston 1965, the
  // rule behind the European rolled-section tables; ~1-3 % of the tables).
  const double tw = s.tw, tf = s.tf;
  const double alpha = -0.042 + 0.2204 * tw / tf + 0.1355 * r / tf - 0.0865 * r * tw / (tf * tf) -
                       0.0725 * tw * tw / (tf * tf);
  const double D = ((tf + r) * (tf + r) + tw * (r + tw / 4.0)) / (2 * r + tf);
  p.J = 2.0 / 3.0 * s.b * std::pow(tf, 3) + (s.h - 2 * tf) * std::pow(tw, 3) / 3.0 +
        2 * alpha * std::pow(D, 4) - 0.420 * std::pow(tf, 4);
  // Shear area along z, EN 1993-1-1 §6.2.6(3)a for rolled I/H sections.
  p.Avz = std::max(p.A - 2 * s.b * tf + (tw + 2 * r) * tf, hw * tw);
  return p;
}

SectionProps props_box(const ShapeBox& s) {
  check_pos(s.h, "h");
  check_pos(s.b, "b");
  check_pos(s.tw, "tw");
  check_pos(s.tf, "tf");
  if (2 * s.tf >= s.h || 2 * s.tw >= s.b) fail("section: box walls overlap");
  SectionProps p;
  double bi = s.b - 2 * s.tw, hi = s.h - 2 * s.tf;
  p.A = s.b * s.h - bi * hi;
  p.Iy = (s.b * std::pow(s.h, 3) - bi * std::pow(hi, 3)) / 12.0;
  p.Iz = (s.h * std::pow(s.b, 3) - hi * std::pow(bi, 3)) / 12.0;
  // Closed thin-wall torsion (Bredt): J = 4 A_m^2 / (integral ds/t)
  double bm = s.b - s.tw, hm = s.h - s.tf;
  double Am = bm * hm;
  p.J = 4.0 * Am * Am / (2.0 * bm / s.tf + 2.0 * hm / s.tw);
  p.Avz = 2 * hi * s.tw;
  p.Avy = 2 * bi * s.tf;
  p.Wely = p.Iy / (s.h / 2);
  p.Welz = p.Iz / (s.b / 2);
  p.Wply = s.b * s.tf * (s.h - s.tf) + s.tw * hi * hi / 2.0;
  p.Wplz = s.h * s.tw * (s.b - s.tw) + s.tf * bi * bi / 2.0;
  return p;
}

SectionProps props_t(const ShapeT& s) {
  check_pos(s.h, "h");
  check_pos(s.b, "b");
  check_pos(s.tw, "tw");
  check_pos(s.tf, "tf");
  if (s.tf >= s.h) fail("section: T flange thicker than height");
  SectionProps p;
  double hw = s.h - s.tf;
  double Af = s.b * s.tf, Aw = hw * s.tw;
  p.A = Af + Aw;
  // z measured from the web bottom, flange on top: centroid height
  double zf = hw + s.tf / 2, zw = hw / 2;
  double zc = (Af * zf + Aw * zw) / p.A;
  p.Iy = s.b * std::pow(s.tf, 3) / 12.0 + Af * (zf - zc) * (zf - zc) +
         s.tw * std::pow(hw, 3) / 12.0 + Aw * (zw - zc) * (zw - zc);
  p.Iz = s.tf * std::pow(s.b, 3) / 12.0 + hw * std::pow(s.tw, 3) / 12.0;
  p.J = (s.b * std::pow(s.tf, 3) + hw * std::pow(s.tw, 3)) / 3.0;
  p.Avz = hw * s.tw;
  p.Avy = Af * 5. / 6.;
  p.Wely = p.Iy / std::max(zc, s.h - zc);
  p.Welz = p.Iz / (s.b / 2);
  return p;
}

SectionProps props_polygon(const ShapePolygon& s) {
  if (s.outer.size() < 3) fail("section: polygon needs >= 3 vertices");
  // Green's theorem integrals over a ring (y,z), CCW positive.
  auto ring_integrals = [](const std::vector<Vec3>& ring, double sign, double& A, double& Sy,
                           double& Sz, double& Iyy, double& Izz, double& Iyz) {
    std::size_t n = ring.size();
    for (std::size_t i = 0; i < n; ++i) {
      double y0 = ring[i].x(), z0 = ring[i].y();
      double y1 = ring[(i + 1) % n].x(), z1 = ring[(i + 1) % n].y();
      double cr = y0 * z1 - y1 * z0;
      A += sign * cr / 2.0;
      Sy += sign * cr * (z0 + z1) / 6.0;                              // integral z dA
      Sz += sign * cr * (y0 + y1) / 6.0;                              // integral y dA
      Iyy += sign * cr * (z0 * z0 + z0 * z1 + z1 * z1) / 12.0;        // integral z^2 dA
      Izz += sign * cr * (y0 * y0 + y0 * y1 + y1 * y1) / 12.0;        // integral y^2 dA
      Iyz += sign * cr * (y0 * z1 + 2 * y0 * z0 + 2 * y1 * z1 + y1 * z0) / 24.0;
    }
  };
  auto signed_area = [](const std::vector<Vec3>& ring) {
    double a = 0;
    std::size_t n = ring.size();
    for (std::size_t i = 0; i < n; ++i)
      a += ring[i].x() * ring[(i + 1) % n].y() - ring[(i + 1) % n].x() * ring[i].y();
    return a / 2.0;
  };

  double A = 0, Sy = 0, Sz = 0, Iyy = 0, Izz = 0, Iyz = 0;
  std::vector<Vec3> outer = s.outer;
  if (signed_area(outer) < 0) std::reverse(outer.begin(), outer.end());
  ring_integrals(outer, +1.0, A, Sy, Sz, Iyy, Izz, Iyz);
  for (auto hole : s.holes) {
    if (signed_area(hole) < 0) std::reverse(hole.begin(), hole.end());
    ring_integrals(hole, -1.0, A, Sy, Sz, Iyy, Izz, Iyz);
  }
  if (!(A > 0) || !std::isfinite(A)) fail("section: degenerate polygon (non-positive area)");

  double yc = Sz / A, zc = Sy / A;
  SectionProps p;
  p.A = A;
  p.Iy = Iyy - A * zc * zc;
  p.Iz = Izz - A * yc * yc;
  double Pyz = Iyz - A * yc * zc;
  // The frame element assumes principal axes; a strongly non-principal polygon
  // is a modeling error we surface instead of silently mis-bending.
  if (std::fabs(Pyz) > 0.05 * std::max(p.Iy, p.Iz))
    fail("section: polygon axes are not principal (|Iyz| > 5% of max(Iy,Iz)); rotate the input");
  // Extreme fibers relative to centroid.
  double ymax = 0, zmax = 0;
  for (const auto& v : outer) {
    ymax = std::max(ymax, std::fabs(v.x() - yc));
    zmax = std::max(zmax, std::fabs(v.y() - zc));
  }
  p.Wely = zmax > 0 ? p.Iy / zmax : 0;
  p.Welz = ymax > 0 ? p.Iz / ymax : 0;
  // Compact-section heuristics for torsion/shear (documented approximations —
  // an exact warping FE is out of scope for the polygon path).
  p.J = std::min(p.Iy + p.Iz, 4.0 * A * A * A * A / (40.0 * (p.Iy + p.Iz)));  // Saint-Venant bound
  p.Avy = p.Avz = A * 5.0 / 6.0;
  return p;
}

namespace {
std::atomic<int> g_shear_areas{0};
}
void set_default_shear_areas(ShearAreas s) { g_shear_areas = s == ShearAreas::Energy ? 0 : 1; }
ShearAreas default_shear_areas() { return g_shear_areas == 0 ? ShearAreas::Energy : ShearAreas::ShapeFactor; }

void Section::finalize() {
  struct Visitor {
    Section& s;
    void operator()(const ShapeGeneric&) {
      require(s.props.A > 0 && s.props.Iy >= 0 && s.props.Iz >= 0,
              "section '" + s.name + "': generic section needs A > 0, Iy/Iz >= 0");
    }
    void operator()(const ShapeRect& sh) {
      s.props = props_rect(sh);
      s.stress_points = {{"+y+z", sh.b / 2, sh.h / 2, 0, 0},
                         {"-y+z", -sh.b / 2, sh.h / 2, 0, 0},
                         {"+y-z", sh.b / 2, -sh.h / 2, 0, 0},
                         {"-y-z", -sh.b / 2, -sh.h / 2, 0, 0},
                         {"center", 0, 0, 1, 1}};
    }
    void operator()(const ShapeCircle& sh) {
      double r = sh.d / 2;
      s.props = props_circle(sh);
      s.stress_points = {{"+y", r, 0, 0, 1}, {"-y", -r, 0, 0, 1}, {"+z", 0, r, 1, 0},
                         {"-z", 0, -r, 1, 0}, {"center", 0, 0, 1, 1}};
    }
    void operator()(const ShapeTube& sh) {
      double r = sh.d / 2;
      s.props = props_tube(sh);
      s.stress_points = {{"+y", r, 0, 0, 1}, {"-y", -r, 0, 0, 1}, {"+z", 0, r, 1, 0},
                         {"-z", 0, -r, 1, 0}};
    }
    void operator()(const ShapeI& sh) {
      s.props = props_i(sh);
      s.stress_points = {{"top+y", sh.b / 2, sh.h / 2, 0, 0},
                         {"top-y", -sh.b / 2, sh.h / 2, 0, 0},
                         {"bot+y", sh.b / 2, -sh.h / 2, 0, 0},
                         {"bot-y", -sh.b / 2, -sh.h / 2, 0, 0},
                         {"web", 0, 0, 0, 1}};
    }
    void operator()(const ShapeBox& sh) {
      s.props = props_box(sh);
      s.stress_points = {{"top+y", sh.b / 2, sh.h / 2, 0, 0},
                         {"top-y", -sh.b / 2, sh.h / 2, 0, 0},
                         {"bot+y", sh.b / 2, -sh.h / 2, 0, 0},
                         {"bot-y", -sh.b / 2, -sh.h / 2, 0, 0},
                         {"web", 0, 0, 0, 1}};
    }
    void operator()(const ShapeT& sh) {
      s.props = props_t(sh);
      double hw = sh.h - sh.tf;
      double Af = sh.b * sh.tf, Aw = hw * sh.tw;
      double zc = (Af * (hw + sh.tf / 2) + Aw * hw / 2) / (Af + Aw);
      s.stress_points = {{"flange+y", sh.b / 2, sh.h - zc, 0, 0},
                         {"flange-y", -sh.b / 2, sh.h - zc, 0, 0},
                         {"stem", 0, -zc, 0, 1}};
    }
    void operator()(const ShapePolygon& sh) {
      s.props = props_polygon(sh);
      // Recenter stress points on the centroid: use extreme vertices.
      double A = 0, Sy = 0, Sz = 0;
      std::size_t n = sh.outer.size();
      for (std::size_t i = 0; i < n; ++i) {
        double y0 = sh.outer[i].x(), z0 = sh.outer[i].y();
        double y1 = sh.outer[(i + 1) % n].x(), z1 = sh.outer[(i + 1) % n].y();
        double cr = y0 * z1 - y1 * z0;
        A += cr / 2;
        Sz += cr * (y0 + y1) / 6.0;
        Sy += cr * (z0 + z1) / 6.0;
      }
      double sgn = A < 0 ? -1.0 : 1.0;
      A *= sgn; Sy *= sgn; Sz *= sgn;
      double yc = Sz / A, zc = Sy / A;
      s.stress_points.clear();
      for (std::size_t i = 0; i < sh.outer.size(); ++i)
        s.stress_points.push_back(
            {"v" + std::to_string(i), sh.outer[i].x() - yc, sh.outer[i].y() - zc, 0, 0});
    }
  };
  std::visit(Visitor{*this}, shape);
  require(props.A > 0, "section '" + name + "': area must be positive");
  if (!std::holds_alternative<ShapeGeneric>(shape)) {
    avy_shape = props.Avy;
    avz_shape = props.Avz;
    shear_nu = std::numeric_limits<double>::quiet_NaN();
  }
}

}  // namespace keta
