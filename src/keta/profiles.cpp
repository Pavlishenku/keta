// Hot-rolled profile outlines — see profiles.h.
#include "keta/profiles.h"

#include <cctype>
#include <cmath>

#include "keta/catalog.h"

namespace keta {

std::vector<Vec3> filleted_polygon(const std::vector<Vec3>& c, const std::vector<double>& r, int seg) {
  require(c.size() >= 3 && r.size() == c.size() && seg >= 1, "filleted polygon: bad input");
  const std::size_t n = c.size();
  std::vector<Vec3> out;
  for (std::size_t k = 0; k < n; ++k) {
    const Vec3 p = c[k], a = c[(k + n - 1) % n], b = c[(k + 1) % n];
    if (r[k] <= 0) {
      out.push_back(p);
      continue;
    }
    const Vec3 u = (a - p).normalized(), v = (b - p).normalized();  // edge directions away from p
    const double th = std::acos(std::max(-1.0, std::min(1.0, u.dot(v))));  // interior angle at p
    const double t = r[k] / std::tan(th / 2);  // tangent length
    require(t < (a - p).norm() && t < (b - p).norm(), "filleted polygon: fillet longer than its edge");
    const Vec3 pa = p + t * u, pb = p + t * v;
    const Vec3 bis = (u + v).normalized();
    const Vec3 cen = p + bis * (r[k] / std::sin(th / 2));
    double a0 = std::atan2(pa.y() - cen.y(), pa.x() - cen.x());
    double a1 = std::atan2(pb.y() - cen.y(), pb.x() - cen.x());
    double da = a1 - a0;
    while (da > kPi) da -= 2 * kPi;
    while (da < -kPi) da += 2 * kPi;
    for (int i = 0; i <= seg; ++i) {
      const double ang = a0 + da * i / seg;
      out.emplace_back(cen.x() + r[k] * std::cos(ang), cen.y() + r[k] * std::sin(ang), 0.0);
    }
  }
  return out;
}

PolygonInertia polygon_inertia(const std::vector<Vec3>& ring) {
  double A = 0, Sy = 0, Sz = 0, Iyy = 0, Izz = 0, Iyz = 0;
  const std::size_t n = ring.size();
  for (std::size_t i = 0; i < n; ++i) {
    const double y0 = ring[i].x(), z0 = ring[i].y(), y1 = ring[(i + 1) % n].x(), z1 = ring[(i + 1) % n].y();
    const double cr = y0 * z1 - y1 * z0;
    A += cr / 2;
    Sy += cr * (z0 + z1) / 6;
    Sz += cr * (y0 + y1) / 6;
    Iyy += cr * (z0 * z0 + z0 * z1 + z1 * z1) / 12;
    Izz += cr * (y0 * y0 + y0 * y1 + y1 * y1) / 12;
    Iyz += cr * (y0 * z1 + 2 * y0 * z0 + 2 * y1 * z1 + y1 * z0) / 24;
  }
  require(A != 0, "polygon inertia: zero area");
  const double s = A < 0 ? -1.0 : 1.0;  // orientation
  PolygonInertia p;
  p.A = s * A;
  p.yc = Sz / A;
  p.zc = Sy / A;
  p.Iy = s * (Iyy - A * p.zc * p.zc);
  p.Iz = s * (Izz - A * p.yc * p.yc);
  p.Iyz = s * (Iyz - A * p.yc * p.zc);
  const double m = 0.5 * (p.Iy + p.Iz), d = std::hypot(0.5 * (p.Iy - p.Iz), p.Iyz);
  p.I1 = m + d;
  p.I2 = m - d;
  p.alpha = 0.5 * std::atan2(-2 * p.Iyz, p.Iy - p.Iz);
  return p;
}

std::vector<Vec3> outline_angle(double h, double t, double r1, double r2, int seg) {
  const std::vector<Vec3> c = {{0, 0, 0}, {h, 0, 0}, {h, t, 0}, {t, t, 0}, {t, h, 0}, {0, h, 0}};
  return filleted_polygon(c, {0, 0, r2, r1, r2, 0}, seg);
}

std::vector<Vec3> outline_inp(double h, double b, double tw, double tf, double r1, double r2, double slope,
                              int seg) {
  // flange inner face: thickness tf at y = b/4 from the toe, rising toward the web
  const double yq = b / 2 - b / 4;  // abscissa (from the web axis) of the tf measurement
  auto tin = [&](double y) { return tf + slope * (yq - y); };  // flange thickness at |y|
  const double ytoe = b / 2, yweb = tw / 2;
  // bottom half, CCW from the bottom-left toe; then mirrored for the top
  const std::vector<Vec3> c = {{-ytoe, 0, 0},
                               {ytoe, 0, 0},
                               {ytoe, tin(ytoe), 0},
                               {yweb, tin(yweb), 0},
                               {yweb, h - tin(yweb), 0},
                               {ytoe, h - tin(ytoe), 0},
                               {ytoe, h, 0},
                               {-ytoe, h, 0},
                               {-ytoe, h - tin(ytoe), 0},
                               {-yweb, h - tin(yweb), 0},
                               {-yweb, tin(yweb), 0},
                               {-ytoe, tin(ytoe), 0}};
  return filleted_polygon(c, {0, 0, r2, r1, r1, r2, 0, 0, r2, r1, r1, r2}, seg);
}

std::vector<Vec3> outline_upe(double h, double b, double tw, double tf, double r, int seg) {
  const std::vector<Vec3> c = {{0, 0, 0}, {b, 0, 0}, {b, tf, 0}, {tw, tf, 0},
                               {tw, h - tf, 0}, {b, h - tf, 0}, {b, h, 0}, {0, h, 0}};
  return filleted_polygon(c, {0, 0, 0, r, r, 0, 0, 0}, seg);
}

std::optional<std::vector<Vec3>> profile_outline(const std::string& designation) {
  const std::string k = normalize_designation(designation);
  // EN 10056-1: L 100 x 100 x 10, r1 = 12, r2 = 6 mm
  if (k == "L100X10" || k == "L100X100X10") return outline_angle(0.100, 0.010, 0.012, 0.006);
  // DIN 1025-1: I 100 (INP 100), h 100, b 50, tw 4.5, tf 6.8, r1 4.5, r2 2.7 mm, slope 14 %
  if (k == "I100" || k == "INP100") return outline_inp(0.100, 0.050, 0.0045, 0.0068, 0.0045, 0.0027);
  // DIN 1026-2: UPE 100, h 100, b 55, tw 4.5, tf 7.5, r 10 mm
  if (k == "UPE100") return outline_upe(0.100, 0.055, 0.0045, 0.0075, 0.010);
  return std::nullopt;
}

}  // namespace keta
