#include "keta/io.h"

#include <fstream>
#include <sstream>

namespace keta {

namespace {
json::Value vec3(const Vec3& v) { return json::Array{v.x(), v.y(), v.z()}; }
Vec3 to_vec3(const json::Value& a) {
  const auto& arr = a.as_array();
  require(arr.size() == 3, "json: expected a 3-vector");
  return Vec3(arr[0].as_number(), arr[1].as_number(), arr[2].as_number());
}
}  // namespace

json::Value model_to_json(const Model& m) {
  json::Value root = json::Value::object();
  {
    json::Value arr = json::Value::array();
    for (const Node& n : m.nodes) {
      json::Value o = json::Value::object();
      o.set("id", n.id);
      o.set("x", n.x.x());
      o.set("y", n.x.y());
      o.set("z", n.x.z());
      arr.push_back(std::move(o));
    }
    root.set("nodes", std::move(arr));
  }
  {
    json::Value arr = json::Value::array();
    for (const Material& mt : m.materials) {
      json::Value o = json::Value::object();
      o.set("id", mt.id);
      o.set("name", mt.name);
      o.set("E", mt.E);
      o.set("nu", mt.nu);
      o.set("rho", mt.rho);
      o.set("alpha", mt.alpha);
      if (mt.fy > 0) o.set("fy", mt.fy);
      if (mt.G_given > 0) o.set("G", mt.G_given);
      if (mt.concrete) {
        json::Value c = json::Value::object();
        c.set("fck", mt.concrete->fck);
        c.set("rh", mt.concrete->rh);
        c.set("h0", mt.concrete->h0);
        c.set("cement_s", mt.concrete->cement_s);
        o.set("concrete", std::move(c));
      }
      arr.push_back(std::move(o));
    }
    root.set("materials", std::move(arr));
  }
  {
    json::Value arr = json::Value::array();
    for (const Section& s : m.sections) {
      json::Value o = json::Value::object();
      o.set("id", s.id);
      o.set("name", s.name);
      struct V {
        json::Value& o;
        void operator()(const ShapeGeneric&) { o.set("type", "generic"); }
        void operator()(const ShapeRect& r) {
          o.set("type", "rect");
          o.set("b", r.b);
          o.set("h", r.h);
        }
        void operator()(const ShapeCircle& c) {
          o.set("type", "circle");
          o.set("d", c.d);
        }
        void operator()(const ShapeTube& t) {
          o.set("type", "tube");
          o.set("d", t.d);
          o.set("t", t.t);
        }
        void operator()(const ShapeI& i) {
          o.set("type", "I");
          o.set("h", i.h);
          o.set("b", i.b);
          o.set("tw", i.tw);
          o.set("tf", i.tf);
          if (i.r != 0) o.set("r", i.r);
        }
        void operator()(const ShapeBox& b) {
          o.set("type", "box");
          o.set("h", b.h);
          o.set("b", b.b);
          o.set("tw", b.tw);
          o.set("tf", b.tf);
        }
        void operator()(const ShapeT& t) {
          o.set("type", "T");
          o.set("h", t.h);
          o.set("b", t.b);
          o.set("tw", t.tw);
          o.set("tf", t.tf);
        }
        void operator()(const ShapePolygon& p) {
          o.set("type", "polygon");
          json::Value outer = json::Value::array();
          for (const Vec3& v : p.outer) outer.push_back(json::Array{v.x(), v.y()});
          o.set("outer", std::move(outer));
        }
      };
      std::visit(V{o}, s.shape);
      json::Value props = json::Value::object();
      props.set("A", s.props.A);
      props.set("Avy", s.props.Avy);
      props.set("Avz", s.props.Avz);
      props.set("Iy", s.props.Iy);
      props.set("Iz", s.props.Iz);
      props.set("J", s.props.J);
      if (s.props.Iw != 0) props.set("Iw", s.props.Iw);
      if (s.props.beta_x != 0) props.set("beta_x", s.props.beta_x);
      o.set("props", std::move(props));
      arr.push_back(std::move(o));
    }
    root.set("sections", std::move(arr));
  }
  {
    json::Value arr = json::Value::array();
    for (const Bar& b : m.bars) {
      json::Value o = json::Value::object();
      o.set("id", b.id);
      o.set("ni", b.ni);
      o.set("nj", b.nj);
      o.set("material", b.material);
      o.set("section", b.section);
      o.set("kind", b.kind == BarKind::Frame ? "frame" : (b.kind == BarKind::Truss ? "truss" : "cable"));
      if (b.roll != 0) o.set("roll", b.roll);
      bool any_rel = false;
      for (bool r : b.releases) any_rel |= r;
      if (any_rel) {
        json::Value rel = json::Value::array();
        for (bool r : b.releases) rel.push_back(r);
        o.set("releases", std::move(rel));
      }
      if (b.has_ecc()) {
        o.set("ecc_i", vec3(b.ecc[0]));
        o.set("ecc_j", vec3(b.ecc[1]));
      }
      if (b.soil.any()) {
        json::Value s = json::Value::object();
        s.set("kfx", b.soil.kfx);
        s.set("kfy", b.soil.kfy);
        s.set("kfz", b.soil.kfz);
        s.set("kmx", b.soil.kmx);
        o.set("soil", std::move(s));
      }
      if (b.cable) {
        json::Value c = json::Value::object();
        c.set("pretension", b.cable->pretension);
        c.set("unstressed_length", b.cable->unstressed_length);
        // appended keys, written only when set (older files stay byte-identical)
        if (b.cable->ea > 0) c.set("ea", b.cable->ea);
        if (b.cable->w) c.set("w", *b.cable->w);
        if (b.cable->end_tension > 0) c.set("end_tension", b.cable->end_tension);
        o.set("cable", std::move(c));
      }
      arr.push_back(std::move(o));
    }
    root.set("bars", std::move(arr));
  }
  // Flat shells: the descriptive slab records and the lowered facets. Both
  // blocks are omitted entirely for a bar-only model (byte-identical output).
  if (!m.slabs.empty()) {
    json::Value arr = json::Value::array();
    for (const Slab& s : m.slabs) {
      json::Value o = json::Value::object();
      o.set("id", s.id);
      o.set("t", s.t);
      o.set("material", s.material);
      o.set("mesh_size", s.mesh_size);
      json::Value ring = json::Value::array();
      for (const Vec3& p : s.outline) ring.push_back(vec3(p));
      o.set("outline", std::move(ring));
      if (!s.holes.empty()) {
        json::Value hs = json::Value::array();
        for (const auto& h : s.holes) {
          json::Value r = json::Value::array();
          for (const Vec3& p : h) r.push_back(vec3(p));
          hs.push_back(std::move(r));
        }
        o.set("holes", std::move(hs));
      }
      if (s.local_x) o.set("local_x", vec3(*s.local_x));
      arr.push_back(std::move(o));
    }
    root.set("slabs", std::move(arr));
  }
  if (!m.shells.empty()) {
    json::Value arr = json::Value::array();
    for (const ShellTri& sh : m.shells) {
      json::Value o = json::Value::object();
      o.set("id", sh.id);
      o.set("n", json::Array{sh.n[0], sh.n[1], sh.n[2]});
      o.set("t", sh.t);
      o.set("material", sh.material);
      o.set("slab", sh.slab);
      if (!sh.enrich.plain()) {
        const ShellEnrich& e = sh.enrich;
        if (e.tn[0] != 0 || e.tn[1] != 0 || e.tn[2] != 0) o.set("t_corners", json::Array{e.tn[0], e.tn[1], e.tn[2]});
        if (e.ecc != 0) o.set("offset", e.ecc);
        if (e.kbed != 0) o.set("bed", e.kbed);
        if (e.ortho) {
          const ShellOrtho& r = *e.ortho;
          json::Value ob = json::Value::object();
          ob.set("Ex", r.Ex);
          ob.set("Ey", r.Ey);
          ob.set("nu_xy", r.nu_xy);
          ob.set("Gxy", r.Gxy);
          ob.set("Gxz", r.Gxz);
          ob.set("Gyz", r.Gyz);
          ob.set("angle", r.angle);
          o.set("orthotropic", std::move(ob));
        }
        if (e.laminate) {
          json::Value lb = json::Value::object();
          json::Value ps = json::Value::array();
          for (const ShellPly& p : e.laminate->plies) {
            json::Value po = json::Value::object();
            po.set("t", p.t);
            po.set("E1", p.E1);
            po.set("E2", p.E2);
            po.set("nu12", p.nu12);
            po.set("G12", p.G12);
            po.set("G13", p.G13);
            po.set("G23", p.G23);
            po.set("angle", p.angle);
            ps.push_back(std::move(po));
          }
          lb.set("plies", std::move(ps));
          if (!e.laminate->centred) lb.set("z_bottom", e.laminate->z_bottom);
          o.set("laminate", std::move(lb));
        }
      }
      arr.push_back(std::move(o));
    }
    root.set("shells", std::move(arr));
  }
  if (!m.continua.empty()) {
    json::Value arr = json::Value::array();
    for (const Continuum& c : m.continua) {
      json::Value o = json::Value::object();
      o.set("id", c.id);
      const char* kinds[] = {"plane_stress", "plane_strain", "axisymmetric", "solid"};
      o.set("kind", kinds[static_cast<int>(c.kind)]);
      json::Value nn = json::Value::array();
      for (Id n : c.n) nn.push_back(n);
      o.set("n", std::move(nn));
      o.set("material", c.material);
      o.set("t", c.t);
      arr.push_back(std::move(o));
    }
    root.set("continua", std::move(arr));
  }
  {
    json::Value arr = json::Value::array();
    for (const Support& s : m.supports) {
      json::Value o = json::Value::object();
      o.set("node", s.node);
      json::Value f = json::Value::array(), k = json::Value::array(), os = json::Value::array();
      for (int c = 0; c < 6; ++c) {
        f.push_back(s.fixed[static_cast<std::size_t>(c)]);
        k.push_back(s.spring[static_cast<std::size_t>(c)]);
        os.push_back(s.oneside[static_cast<std::size_t>(c)]);
      }
      o.set("fixed", std::move(f));
      o.set("spring", std::move(k));
      o.set("oneside", std::move(os));
      bool any_cap = false;
      for (double cp : s.cap) any_cap |= cp != 0;
      if (any_cap) {
        json::Value cp = json::Value::array();
        for (double x : s.cap) cp.push_back(x);
        o.set("cap", std::move(cp));
      }
      arr.push_back(std::move(o));
    }
    root.set("supports", std::move(arr));
  }
  {
    json::Value arr = json::Value::array();
    for (const PointMass& pm : m.point_masses) {
      json::Value o = json::Value::object();
      o.set("node", pm.node);
      o.set("mass", pm.mass);
      o.set("inertia", vec3(pm.inertia));
      arr.push_back(std::move(o));
    }
    if (!m.point_masses.empty()) root.set("point_masses", std::move(arr));
  }
  if (!m.tendons.empty()) {
    // bar-chain prestressing tendons: every field, so a model read back
    // tensions the same tendons in a phasing
    auto arr4 = [](const std::array<double, 4>& a) {
      json::Value v = json::Value::array();
      for (double x : a) v.push_back(x);
      return v;
    };
    json::Value arr = json::Value::array();
    for (const Tendon& t : m.tendons) {
      json::Value o = json::Value::object();
      o.set("id", t.id);
      o.set("name", t.name);
      json::Value bars = json::Value::array();
      for (Id b : t.bars) bars.push_back(b);
      o.set("bars", std::move(bars));
      json::Value prof = json::Value::array();
      for (const PolySegment& g : t.profile) {
        json::Value q = json::Value::object();
        q.set("s0", g.s0);
        q.set("s1", g.s1);
        q.set("cy", arr4(g.cy));
        q.set("cz", arr4(g.cz));
        if (g.exact != 0) {
          q.set("exact", g.exact);
          json::Value pp = json::Value::array();
          for (const Vec3& x : g.p) pp.push_back(vec3(x));
          q.set("p", std::move(pp));
          q.set("C", vec3(g.C));
          q.set("U", vec3(g.U));
          q.set("D", vec3(g.D));
          q.set("R", g.R);
          q.set("th", g.th);
        }
        prof.push_back(std::move(q));
      }
      o.set("profile", std::move(prof));
      o.set("jack_force", t.jack_force);
      o.set("jack_both_ends", t.jack_both_ends);
      const TendonLosses& l = t.losses;
      json::Value lo = json::Value::object();
      lo.set("mu", l.mu);
      lo.set("k", l.k);
      lo.set("slip", l.slip);
      lo.set("Ep", l.Ep);
      lo.set("Ap", l.Ap);
      lo.set("rho1000", l.rho1000);
      lo.set("fpk", l.fpk);
      lo.set("relaxation_class", l.relaxation_class);
      lo.set("relaxation_bpel", l.relaxation_bpel);
      lo.set("relax_teb_days", l.relax_teb_days);
      lo.set("thermal_loss", l.thermal_loss);
      lo.set("phi_abs", l.phi_abs);
      lo.set("non_simul", l.non_simul);
      o.set("losses", std::move(lo));
      json::Value fp = json::Value::array();
      for (const auto& [sx, P] : t.force_profile) {
        json::Value pr = json::Value::array();
        pr.push_back(sx);
        pr.push_back(P);
        fp.push_back(std::move(pr));
      }
      o.set("force_profile", std::move(fp));
      json::Value fl = json::Value::array();
      for (bool f : t.flips) fl.push_back(f);
      o.set("flips", std::move(fl));
      o.set("case_id", t.case_id);
      o.set("s_begin", t.s_begin);
      o.set("s_end", t.s_end);
      o.set("exact_angle", t.exact_angle);
      o.set("jack_from_end", t.jack_from_end);
      o.set("release_end_first", t.release_end_first);
      json::Value fr = json::Value::array();
      for (const auto& [a, b] : t.free_runs) {
        json::Value pr = json::Value::array();
        pr.push_back(a);
        pr.push_back(b);
        fr.push_back(std::move(pr));
      }
      o.set("free_runs", std::move(fr));
      arr.push_back(std::move(o));
    }
    root.set("tendons", std::move(arr));
  }
  if (!m.shell_tendons.empty()) {
    json::Value arr = json::Value::array();
    for (const ShellTendon& t : m.shell_tendons) {
      json::Value o = json::Value::object();
      o.set("id", t.id);
      json::Value pts = json::Value::array();
      for (const Vec3& p : t.pts) pts.push_back(vec3(p));
      o.set("points", std::move(pts));
      o.set("P0", t.P0);
      o.set("mu", t.mu);
      o.set("k", t.k);
      o.set("jack", t.jack == ShellTendon::Jack::Start ? "start" : t.jack == ShellTendon::Jack::End ? "end" : "both");
      o.set("slab", t.slab);
      if (t.Ap > 0) {
        o.set("Ap", t.Ap);
        o.set("Ep", t.Ep);
        o.set("fpk", t.fpk);
        o.set("rho1000", t.rho1000);
        o.set("relaxation_class", t.relaxation_class);
      }
      arr.push_back(std::move(o));
    }
    root.set("shell_tendons", std::move(arr));
  }
  {
    json::Value arr = json::Value::array();
    for (const LoadCase& lc : m.cases) {
      json::Value o = json::Value::object();
      o.set("id", lc.id);
      o.set("name", lc.name);
      json::Value loads = json::Value::array();
      for (const Load& l : lc.loads) {
        json::Value lo = json::Value::object();
        const char* frame = l.frame == LoadFrame::Local
                                ? "local"
                                : (l.frame == LoadFrame::GlobalProjected ? "projected" : "global");
        switch (l.kind) {
          case Load::Kind::NodalForce:
            lo.set("type", "nodal");
            lo.set("node", l.node);
            lo.set("F", vec3(l.F));
            lo.set("M", vec3(l.M));
            break;
          case Load::Kind::BarUdl:
            lo.set("type", "udl");
            lo.set("bar", l.bar);
            lo.set("q", vec3(l.F));
            lo.set("frame", frame);
            lo.set("a", l.a);
            lo.set("b", l.b);
            break;
          case Load::Kind::BarLinear:
            lo.set("type", "linear");
            lo.set("bar", l.bar);
            lo.set("qa", vec3(l.F));
            lo.set("qb", vec3(l.F2));
            lo.set("frame", frame);
            lo.set("a", l.a);
            lo.set("b", l.b);
            break;
          case Load::Kind::BarParabolic:
            lo.set("type", "parabolic");
            lo.set("bar", l.bar);
            lo.set("qa", vec3(l.F));
            lo.set("qmid", vec3(l.Fmid));
            lo.set("qb", vec3(l.F2));
            lo.set("frame", frame);
            lo.set("a", l.a);
            lo.set("b", l.b);
            break;
          case Load::Kind::BarPoint:
            lo.set("type", "point");
            lo.set("bar", l.bar);
            lo.set("a", l.a);
            lo.set("F", vec3(l.F));
            lo.set("M", vec3(l.M));
            lo.set("frame", frame);
            break;
          case Load::Kind::BarThermal:
            lo.set("type", "thermal");
            lo.set("bar", l.bar);
            lo.set("dt", l.dt);
            lo.set("gy", l.gy);
            lo.set("gz", l.gz);
            break;
          case Load::Kind::SelfWeight:
            lo.set("type", "self_weight");
            lo.set("g", vec3(l.g));
            break;
          case Load::Kind::SupportDisp: {
            lo.set("type", "support_disp");
            lo.set("node", l.node);
            json::Value d = json::Value::array();
            for (int c = 0; c < 6; ++c)
              d.push_back(l.disp[static_cast<std::size_t>(c)]
                              ? json::Value(*l.disp[static_cast<std::size_t>(c)])
                              : json::Value(nullptr));
            lo.set("disp", std::move(d));
            break;
          }
          case Load::Kind::Prestrain:
            lo.set("type", "prestrain");
            lo.set("bar", l.bar);
            lo.set("strain", l.strain);
            lo.set("kappa_y", l.kappa_y);
            lo.set("kappa_z", l.kappa_z);
            if (l.kappa_x != 0) lo.set("kappa_x", l.kappa_x);
            break;
          case Load::Kind::PrestrainPoint:
            lo.set("type", "prestrain_point");
            lo.set("bar", l.bar);
            lo.set("a", l.a);
            lo.set("dux", l.strain);
            lo.set("theta_y", l.kappa_y);
            lo.set("theta_z", l.kappa_z);
            if (l.dvy != 0) lo.set("dvy", l.dvy);
            if (l.dvz != 0) lo.set("dvz", l.dvz);
            if (l.dthx != 0) lo.set("theta_x", l.dthx);
            if (l.frame == LoadFrame::Global) lo.set("frame", "global");
            break;
          case Load::Kind::PrestrainLinear:
            lo.set("type", "prestrain_linear");
            lo.set("bar", l.bar);
            lo.set("strain_i", l.strain);
            lo.set("kappa_y_i", l.kappa_y);
            lo.set("kappa_z_i", l.kappa_z);
            lo.set("strain_j", l.strain2);
            lo.set("kappa_y_j", l.kappa_y2);
            lo.set("kappa_z_j", l.kappa_z2);
            if (l.kappa_x != 0 || l.kappa_x2 != 0) {
              lo.set("kappa_x_i", l.kappa_x);
              lo.set("kappa_x_j", l.kappa_x2);
            }
            break;
          case Load::Kind::PrestrainParabolic:
            lo.set("type", "prestrain_parabolic");
            lo.set("bar", l.bar);
            lo.set("strain_i", l.strain);
            lo.set("kappa_y_i", l.kappa_y);
            lo.set("kappa_z_i", l.kappa_z);
            lo.set("kappa_x_i", l.kappa_x);
            lo.set("strain_m", l.strain_m);
            lo.set("kappa_y_m", l.kappa_y_m);
            lo.set("kappa_z_m", l.kappa_z_m);
            lo.set("kappa_x_m", l.kappa_x_m);
            lo.set("strain_j", l.strain2);
            lo.set("kappa_y_j", l.kappa_y2);
            lo.set("kappa_z_j", l.kappa_z2);
            lo.set("kappa_x_j", l.kappa_x2);
            break;
          case Load::Kind::SurfacePressure:
            lo.set("type", "surface_pressure");
            lo.set("slab", l.slab);
            lo.set("w", vec3(l.F));
            lo.set("frame", l.sframe == SurfaceFrame::LocalNormal
                                ? "local_normal"
                                : (l.sframe == SurfaceFrame::Projected ? "projected" : "global"));
            break;
          case Load::Kind::ShellPrestress:
            lo.set("type", "shell_prestress");
            lo.set("tendon", l.node);
            lo.set("factor", l.strain);
            break;
          case Load::Kind::ShellImposed: {
            lo.set("type", "shell_imposed");
            lo.set("shell", l.node);
            json::Value d = json::Value::array();
            for (const Vec3& v : l.poly)
              for (int c = 0; c < 3; ++c) d.push_back(v(c));
            lo.set("d", std::move(d));
            break;
          }
          case Load::Kind::PolygonPressure: {
            lo.set("type", "polygon_pressure");
            lo.set("slab", l.slab);
            lo.set("q", vec3(l.F));
            lo.set("proj", vec3(l.proj));
            json::Value pts = json::Value::array();
            for (const Vec3& p : l.poly) pts.push_back(vec3(p));
            lo.set("polygon", std::move(pts));
            json::Value hs = json::Value::array();
            for (const auto& h : l.holes) {
              json::Value r = json::Value::array();
              for (const Vec3& p : h) r.push_back(vec3(p));
              hs.push_back(std::move(r));
            }
            lo.set("holes", std::move(hs));
            break;
          }
          case Load::Kind::Hydrostatic:
            lo.set("type", "hydrostatic");
            lo.set("slab", l.slab);
            lo.set("gamma", l.hydro_gamma);
            lo.set("point", vec3(l.hydro_point));
            lo.set("down", vec3(l.g));
            lo.set("side", l.hydro_side);
            break;
          case Load::Kind::SurfaceThermal:
            lo.set("type", "surface_thermal");
            lo.set("slab", l.slab);
            lo.set("dt", l.dt);
            lo.set("grad", l.grad);
            break;
        }
        // a/b as fractions of the bar length (the Python default): without it
        // the reloaded abscissae would be read in metres
        if (l.relative) lo.set("relative", true);
        loads.push_back(std::move(lo));
      }
      o.set("loads", std::move(loads));
      arr.push_back(std::move(o));
    }
    root.set("cases", std::move(arr));
  }
  {
    json::Value arr = json::Value::array();
    for (const Combination& cb : m.combos) {
      json::Value o = json::Value::object();
      o.set("id", cb.id);
      o.set("name", cb.name);
      json::Value terms = json::Value::array();
      for (const auto& [cid, f] : cb.terms) terms.push_back(json::Array{cid, f});
      o.set("terms", std::move(terms));
      arr.push_back(std::move(o));
    }
    if (!m.combos.empty()) root.set("combinations", std::move(arr));
  }
  // Links, linear constraints, skew supports and the added masses: written
  // so that a model read back gives the same results (they were silently
  // dropped before). Tapered and composite bars hold computed kernels
  // (closures), not data: such a model is refused rather than written
  // without them.
  require(m.tapered.empty(),
          "json: a model with tapered bars (set_tapered) cannot be written: their kernel is not serialisable");
  require(m.composite.empty(),
          "json: a model with composite bars cannot be written: their parts are not serialisable");
  auto flat = [](const auto& mat) {
    json::Value a = json::Value::array();
    for (Eigen::Index r = 0; r < mat.rows(); ++r)
      for (Eigen::Index c = 0; c < mat.cols(); ++c) a.push_back(mat(r, c));
    return a;
  };
  auto bools6 = [](const std::array<bool, 6>& b) {
    json::Value a = json::Value::array();
    for (bool x : b) a.push_back(x);
    return a;
  };
  if (!m.rigid_links.empty()) {
    json::Value arr = json::Value::array();
    for (const RigidLink& rl : m.rigid_links) {
      json::Value o = json::Value::object();
      o.set("id", rl.id);
      o.set("master", rl.master);
      json::Value s = json::Value::array();
      for (Id x : rl.slaves) s.push_back(x);
      o.set("slaves", std::move(s));
      o.set("dofs", bools6(rl.dofs));
      arr.push_back(std::move(o));
    }
    root.set("rigid_links", std::move(arr));
  }
  if (!m.links.empty()) {
    json::Value arr = json::Value::array();
    for (const ElasticLink& el : m.links) {
      json::Value o = json::Value::object();
      o.set("id", el.id);
      o.set("ni", el.ni);
      o.set("nj", el.nj);
      o.set("k", flat(el.k));
      if (el.axes) o.set("axes", flat(*el.axes));
      o.set("roll", el.roll);
      o.set("s", el.s);
      if (el.has_law()) {
        json::Value laws = json::Value::array();
        for (const auto& law : el.law) {
          if (!law) {
            laws.push_back(json::Value());
            continue;
          }
          json::Value pts = json::Value::array();
          for (const auto& [d, f] : law->pts) pts.push_back(json::Array{d, f});
          laws.push_back(std::move(pts));
        }
        o.set("laws", std::move(laws));
      }
      arr.push_back(std::move(o));
    }
    root.set("links", std::move(arr));
  }
  if (!m.constraints.empty()) {
    json::Value arr = json::Value::array();
    for (const Constraint& c : m.constraints) {
      json::Value o = json::Value::object();
      o.set("id", c.id);
      json::Value terms = json::Value::array();
      for (const ConstraintTerm& t : c.terms) terms.push_back(json::Array{t.node, t.comp, t.coef});
      o.set("terms", std::move(terms));
      arr.push_back(std::move(o));
    }
    root.set("constraints", std::move(arr));
  }
  if (!m.skew_supports.empty()) {
    json::Value arr = json::Value::array();
    for (const SkewSupport& sk : m.skew_supports) {
      json::Value o = json::Value::object();
      o.set("id", sk.id);
      o.set("node", sk.node);
      o.set("axes", flat(sk.axes));
      o.set("fixed", bools6(sk.fixed));
      arr.push_back(std::move(o));
    }
    root.set("skew_supports", std::move(arr));
  }
  if (!m.nodal_masses.empty()) {
    json::Value arr = json::Value::array();
    for (const NodalMass& nm : m.nodal_masses) {
      json::Value o = json::Value::object();
      o.set("node", nm.node);
      o.set("m", flat(nm.m));
      if (nm.add) o.set("add", true);
      arr.push_back(std::move(o));
    }
    root.set("nodal_masses", std::move(arr));
  }
  if (!m.bar_masses.empty()) {
    static const char* kinds[] = {"auto", "nul", "uni", "lin", "parabol"};
    static const char* forms[] = {"model", "diagonal", "consistent"};
    json::Value arr = json::Value::array();
    for (const BarMass& bm : m.bar_masses) {
      json::Value o = json::Value::object();
      o.set("bar", bm.bar);
      o.set("kind", std::string(kinds[static_cast<int>(bm.kind)]));
      o.set("m1", bm.m1);
      o.set("m2", bm.m2);
      o.set("m3", bm.m3);
      o.set("ecc_or", vec3(bm.ecc_or));
      o.set("ecc_ex", vec3(bm.ecc_ex));
      o.set("form", std::string(forms[static_cast<int>(bm.form)]));
      arr.push_back(std::move(o));
    }
    root.set("bar_masses", std::move(arr));
  }
  if (m.bar_mass_selective) root.set("bar_mass_selective", true);
  return root;
}

Model model_from_json(const json::Value& v) {
  Model m;
  for (const auto& n : v.at("nodes").as_array()) {
    Node node;
    node.id = n.at("id").as_int();
    node.x = Vec3(n.num_or("x", 0), n.num_or("y", 0), n.num_or("z", 0));
    m.nodes.push_back(node);
  }
  for (const auto& o : v.at("materials").as_array()) {
    Material mt;
    mt.id = o.at("id").as_int();
    mt.name = o.str_or("name", "");
    mt.E = o.at("E").as_number();
    mt.nu = o.num_or("nu", 0.3);
    mt.rho = o.num_or("rho", 0.0);
    mt.alpha = o.num_or("alpha", 1.2e-5);
    mt.fy = o.num_or("fy", 0.0);
    mt.G_given = o.num_or("G", 0.0);
    if (const json::Value* c = o.get("concrete")) {
      ConcreteTime ct;
      ct.fck = c->num_or("fck", 35e6);
      ct.rh = c->num_or("rh", 70);
      ct.h0 = c->num_or("h0", 0.3);
      ct.cement_s = c->num_or("cement_s", 0.25);
      mt.concrete = ct;
    }
    m.materials.push_back(mt);
  }
  for (const auto& o : v.at("sections").as_array()) {
    Section s;
    s.id = o.at("id").as_int();
    s.name = o.str_or("name", "");
    std::string type = o.str_or("type", "generic");
    if (type == "rect")
      s.shape = ShapeRect{o.at("b").as_number(), o.at("h").as_number()};
    else if (type == "circle")
      s.shape = ShapeCircle{o.at("d").as_number()};
    else if (type == "tube")
      s.shape = ShapeTube{o.at("d").as_number(), o.at("t").as_number()};
    else if (type == "I")
      s.shape = ShapeI{o.at("h").as_number(), o.at("b").as_number(), o.at("tw").as_number(),
                       o.at("tf").as_number(), o.num_or("r", 0.0)};
    else if (type == "box")
      s.shape = ShapeBox{o.at("h").as_number(), o.at("b").as_number(), o.at("tw").as_number(),
                         o.at("tf").as_number()};
    else if (type == "T")
      s.shape = ShapeT{o.at("h").as_number(), o.at("b").as_number(), o.at("tw").as_number(),
                       o.at("tf").as_number()};
    else if (type == "polygon") {
      ShapePolygon p;
      for (const auto& pt : o.at("outer").as_array())
        p.outer.push_back(Vec3(pt.as_array()[0].as_number(), pt.as_array()[1].as_number(), 0));
      s.shape = p;
    } else {
      s.shape = ShapeGeneric{};
      const json::Value& props = o.at("props");
      s.props.A = props.at("A").as_number();
      s.props.Avy = props.num_or("Avy", 0);
      s.props.Avz = props.num_or("Avz", 0);
      s.props.Iy = props.num_or("Iy", 0);
      s.props.Iz = props.num_or("Iz", 0);
      s.props.J = props.num_or("J", 0);
      s.props.Iw = props.num_or("Iw", 0);
      s.props.beta_x = props.num_or("beta_x", 0);
    }
    s.finalize();
    m.sections.push_back(s);
  }
  for (const auto& o : v.at("bars").as_array()) {
    Bar b;
    b.id = o.at("id").as_int();
    b.ni = o.at("ni").as_int();
    b.nj = o.at("nj").as_int();
    b.material = o.at("material").as_int();
    b.section = o.at("section").as_int();
    std::string kind = o.str_or("kind", "frame");
    b.kind = kind == "truss" ? BarKind::Truss : (kind == "cable" ? BarKind::Cable : BarKind::Frame);
    b.roll = o.num_or("roll", 0);
    if (const json::Value* rel = o.get("releases")) {
      const auto& arr = rel->as_array();
      for (std::size_t i = 0; i < 12 && i < arr.size(); ++i) b.releases[i] = arr[i].as_bool();
    }
    if (o.has("ecc_i")) b.ecc[0] = to_vec3(o.at("ecc_i"));
    if (o.has("ecc_j")) b.ecc[1] = to_vec3(o.at("ecc_j"));
    if (const json::Value* s = o.get("soil")) {
      b.soil.kfx = s->num_or("kfx", 0);
      b.soil.kfy = s->num_or("kfy", 0);
      b.soil.kfz = s->num_or("kfz", 0);
      b.soil.kmx = s->num_or("kmx", 0);
    }
    if (const json::Value* c = o.get("cable")) {
      CableProps cp;
      cp.pretension = c->num_or("pretension", 0);
      cp.unstressed_length = c->num_or("unstressed_length", 0);
      cp.ea = c->num_or("ea", 0);
      if (const json::Value* w = c->get("w"); w && w->is_number()) cp.w = w->as_number();
      cp.end_tension = c->num_or("end_tension", 0);
      b.cable = cp;
    }
    m.bars.push_back(b);
    m.resolve_shear(m.bars.back());
  }
  // Flat shells (absent in a bar-only document).
  if (const json::Value* sl = v.get("slabs"))
    for (const auto& o : sl->as_array()) {
      Slab s;
      s.id = o.at("id").as_int();
      s.t = o.num_or("t", 0);
      s.material = o.at("material").as_int();
      s.mesh_size = o.num_or("mesh_size", 0);
      for (const auto& p : o.at("outline").as_array()) s.outline.push_back(to_vec3(p));
      if (const json::Value* hs = o.get("holes"))
        for (const auto& h : hs->as_array()) {
          std::vector<Vec3> ring;
          for (const auto& p : h.as_array()) ring.push_back(to_vec3(p));
          s.holes.push_back(std::move(ring));
        }
      if (o.has("local_x")) s.local_x = to_vec3(o.at("local_x"));
      m.slabs.push_back(std::move(s));
    }
  if (const json::Value* sh = v.get("shells"))
    for (const auto& o : sh->as_array()) {
      ShellTri t;
      t.id = o.at("id").as_int();
      const auto& nn = o.at("n").as_array();
      require(nn.size() == 3, "json: shell needs exactly 3 node ids");
      for (std::size_t k = 0; k < 3; ++k) t.n[k] = nn[k].as_int();
      t.t = o.num_or("t", 0);
      t.material = o.at("material").as_int();
      t.slab = o.at("slab").as_int();
      if (const json::Value* tc = o.get("t_corners"))
        for (std::size_t k = 0; k < 3; ++k) t.enrich.tn[k] = tc->as_array()[k].as_number();
      t.enrich.ecc = o.num_or("offset", 0);
      t.enrich.kbed = o.num_or("bed", 0);
      if (const json::Value* ob = o.get("orthotropic")) {
        ShellOrtho r;
        r.Ex = ob->num_or("Ex", 0);
        r.Ey = ob->num_or("Ey", 0);
        r.nu_xy = ob->num_or("nu_xy", 0);
        r.Gxy = ob->num_or("Gxy", 0);
        r.Gxz = ob->num_or("Gxz", 0);
        r.Gyz = ob->num_or("Gyz", 0);
        r.angle = ob->num_or("angle", 0);
        t.enrich.ortho = r;
      }
      if (const json::Value* lb = o.get("laminate")) {
        ShellLaminate lm;
        for (const auto& po : lb->at("plies").as_array()) {
          ShellPly p;
          p.t = po.num_or("t", 0);
          p.E1 = po.num_or("E1", 0);
          p.E2 = po.num_or("E2", 0);
          p.nu12 = po.num_or("nu12", 0);
          p.G12 = po.num_or("G12", 0);
          p.G13 = po.num_or("G13", 0);
          p.G23 = po.num_or("G23", 0);
          p.angle = po.num_or("angle", 0);
          lm.plies.push_back(p);
        }
        if (const json::Value* zb = lb->get("z_bottom")) {
          lm.centred = false;
          lm.z_bottom = zb->as_number();
        }
        t.enrich.laminate = lm;
      }
      m.shells.push_back(t);
    }
  if (const json::Value* cs = v.get("continua"))
    for (const auto& o : cs->as_array()) {
      const std::string kd = o.at("kind").as_string();
      Continuum::Kind k = Continuum::Kind::PlaneStress;
      if (kd == "plane_strain") k = Continuum::Kind::PlaneStrain;
      else if (kd == "axisymmetric") k = Continuum::Kind::Axisymmetric;
      else if (kd == "solid") k = Continuum::Kind::Solid;
      else require(kd == "plane_stress", "json: unknown continuum kind " + kd);
      std::vector<Id> nn;
      for (const auto& x : o.at("n").as_array()) nn.push_back(x.as_int());
      const Id id = m.add_continuum(k, nn, o.at("material").as_int(), o.num_or("t", 1.0));
      m.continua.back().id = o.get("id") ? o.at("id").as_int() : id;
    }
  if (const json::Value* sup = v.get("supports"))
    for (const auto& o : sup->as_array()) {
      Support s;
      s.node = o.at("node").as_int();
      if (const json::Value* f = o.get("fixed"))
        for (int c = 0; c < 6; ++c) s.fixed[static_cast<std::size_t>(c)] = f->as_array()[static_cast<std::size_t>(c)].as_bool();
      if (const json::Value* k = o.get("spring"))
        for (int c = 0; c < 6; ++c) s.spring[static_cast<std::size_t>(c)] = k->as_array()[static_cast<std::size_t>(c)].as_number();
      if (const json::Value* os = o.get("oneside"))
        for (int c = 0; c < 6; ++c)
          s.oneside[static_cast<std::size_t>(c)] = static_cast<int>(os->as_array()[static_cast<std::size_t>(c)].as_int());
      if (const json::Value* cp = o.get("cap"))
        for (int c = 0; c < 6; ++c) s.cap[static_cast<std::size_t>(c)] = cp->as_array()[static_cast<std::size_t>(c)].as_number();
      m.supports.push_back(s);
    }
  if (const json::Value* pm = v.get("point_masses"))
    for (const auto& o : pm->as_array()) {
      PointMass p;
      p.node = o.at("node").as_int();
      p.mass = o.at("mass").as_number();
      if (o.has("inertia")) p.inertia = to_vec3(o.at("inertia"));
      m.point_masses.push_back(p);
    }
  if (const json::Value* td = v.get("tendons"))
    for (const auto& o : td->as_array()) {
      Tendon t;
      t.id = o.at("id").as_int();
      t.name = o.str_or("name", "");
      for (const auto& b : o.at("bars").as_array()) t.bars.push_back(b.as_int());
      auto arr4 = [](const json::Value& v) {
        std::array<double, 4> a{};
        const auto& xs = v.as_array();
        for (std::size_t i = 0; i < 4 && i < xs.size(); ++i) a[i] = xs[i].as_number();
        return a;
      };
      for (const auto& q : o.at("profile").as_array()) {
        PolySegment g;
        g.s0 = q.at("s0").as_number();
        g.s1 = q.at("s1").as_number();
        g.cy = arr4(q.at("cy"));
        g.cz = arr4(q.at("cz"));
        if (q.has("exact")) {
          g.exact = q.at("exact").as_int();
          const auto& pp = q.at("p").as_array();
          for (std::size_t i = 0; i < 4 && i < pp.size(); ++i) g.p[i] = to_vec3(pp[i]);
          g.C = to_vec3(q.at("C"));
          g.U = to_vec3(q.at("U"));
          g.D = to_vec3(q.at("D"));
          g.R = q.num_or("R", 0);
          g.th = q.num_or("th", 0);
        }
        t.profile.push_back(g);
      }
      t.jack_force = o.num_or("jack_force", 0);
      t.jack_both_ends = o.bool_or("jack_both_ends", false);
      if (const json::Value* lo = o.get("losses")) {
        TendonLosses& l = t.losses;
        l.mu = lo->num_or("mu", l.mu);
        l.k = lo->num_or("k", l.k);
        l.slip = lo->num_or("slip", l.slip);
        l.Ep = lo->num_or("Ep", l.Ep);
        l.Ap = lo->num_or("Ap", l.Ap);
        l.rho1000 = lo->num_or("rho1000", l.rho1000);
        l.fpk = lo->num_or("fpk", l.fpk);
        l.relaxation_class = static_cast<int>(lo->num_or("relaxation_class", l.relaxation_class));
        l.relaxation_bpel = lo->bool_or("relaxation_bpel", l.relaxation_bpel);
        l.relax_teb_days = lo->num_or("relax_teb_days", l.relax_teb_days);
        l.thermal_loss = lo->num_or("thermal_loss", l.thermal_loss);
        l.phi_abs = lo->num_or("phi_abs", l.phi_abs);
        l.non_simul = static_cast<int>(lo->num_or("non_simul", l.non_simul));
      }
      if (const json::Value* fp = o.get("force_profile"))
        for (const auto& pr : fp->as_array())
          t.force_profile.emplace_back(pr.as_array().at(0).as_number(), pr.as_array().at(1).as_number());
      if (const json::Value* fl = o.get("flips"))
        for (const auto& f : fl->as_array()) t.flips.push_back(f.as_bool());
      t.case_id = static_cast<Id>(o.num_or("case_id", 0));
      t.s_begin = o.num_or("s_begin", 0.0);
      t.s_end = o.num_or("s_end", -1.0);
      t.exact_angle = o.bool_or("exact_angle", false);
      t.jack_from_end = o.bool_or("jack_from_end", false);
      t.release_end_first = o.bool_or("release_end_first", false);
      if (const json::Value* fr = o.get("free_runs"))
        for (const auto& pr : fr->as_array())
          t.free_runs.emplace_back(pr.as_array().at(0).as_number(), pr.as_array().at(1).as_number());
      m.tendons.push_back(t);
    }
  if (const json::Value* st = v.get("shell_tendons"))
    for (const auto& o : st->as_array()) {
      ShellTendon t;
      t.id = o.at("id").as_int();
      for (const auto& p : o.at("points").as_array()) t.pts.push_back(to_vec3(p));
      t.P0 = o.at("P0").as_number();
      t.mu = o.num_or("mu", 0);
      t.k = o.num_or("k", 0);
      const std::string j = o.str_or("jack", "start");
      t.jack = j == "end" ? ShellTendon::Jack::End : j == "both" ? ShellTendon::Jack::Both : ShellTendon::Jack::Start;
      t.slab = o.has("slab") ? o.at("slab").as_int() : 0;
      t.Ap = o.num_or("Ap", 0);
      t.Ep = o.num_or("Ep", 195e9);
      t.fpk = o.num_or("fpk", 0);
      t.rho1000 = o.num_or("rho1000", 2.5);
      t.relaxation_class = static_cast<int>(o.num_or("relaxation_class", 2));
      m.shell_tendons.push_back(t);
    }
  if (const json::Value* cs = v.get("cases"))
    for (const auto& o : cs->as_array()) {
      LoadCase lc;
      lc.id = o.at("id").as_int();
      lc.name = o.str_or("name", "");
      for (const auto& lo : o.at("loads").as_array()) {
        std::string type = lo.at("type").as_string();
        std::string fr = lo.str_or("frame", "global");
        LoadFrame frame = fr == "local" ? LoadFrame::Local
                          : fr == "projected" ? LoadFrame::GlobalProjected
                                              : LoadFrame::Global;
        Load l;
        if (type == "nodal") {
          l = Load::nodal(lo.at("node").as_int(), to_vec3(lo.at("F")),
                          lo.has("M") ? to_vec3(lo.at("M")) : Vec3::Zero());
        } else if (type == "udl") {
          l = Load::udl(lo.at("bar").as_int(), to_vec3(lo.at("q")), frame, lo.num_or("a", 0),
                        lo.num_or("b", -1));
        } else if (type == "linear") {
          l = Load::linear(lo.at("bar").as_int(), to_vec3(lo.at("qa")), to_vec3(lo.at("qb")),
                           frame, lo.num_or("a", 0), lo.num_or("b", -1));
        } else if (type == "parabolic") {
          l = Load::parabolic(lo.at("bar").as_int(), to_vec3(lo.at("qa")), to_vec3(lo.at("qmid")),
                              to_vec3(lo.at("qb")), frame, lo.num_or("a", 0), lo.num_or("b", -1));
        } else if (type == "point") {
          l = Load::point(lo.at("bar").as_int(), lo.at("a").as_number(), to_vec3(lo.at("F")),
                          lo.has("M") ? to_vec3(lo.at("M")) : Vec3::Zero(), frame);
        } else if (type == "thermal") {
          l = Load::thermal(lo.at("bar").as_int(), lo.num_or("dt", 0), lo.num_or("gy", 0),
                            lo.num_or("gz", 0));
        } else if (type == "self_weight") {
          l = Load::self_weight(lo.has("g") ? to_vec3(lo.at("g")) : Vec3(0, 0, -kGravity));
        } else if (type == "support_disp") {
          std::array<std::optional<double>, 6> d{};
          const auto& arr = lo.at("disp").as_array();
          for (int c = 0; c < 6 && c < static_cast<int>(arr.size()); ++c)
            if (!arr[static_cast<std::size_t>(c)].is_null())
              d[static_cast<std::size_t>(c)] = arr[static_cast<std::size_t>(c)].as_number();
          l = Load::support_disp(lo.at("node").as_int(), d);
        } else if (type == "prestrain") {
          l.kind = Load::Kind::Prestrain;
          l.bar = lo.at("bar").as_int();
          l.strain = lo.num_or("strain", 0);
          l.kappa_y = lo.num_or("kappa_y", 0);
          l.kappa_z = lo.num_or("kappa_z", 0);
          l.kappa_x = lo.num_or("kappa_x", 0);
        } else if (type == "prestrain_point") {
          l = Load::prestrain_point(lo.at("bar").as_int(), lo.at("a").as_number(),
                                    lo.num_or("dux", 0), lo.num_or("theta_y", 0),
                                    lo.num_or("theta_z", 0));
          l.dvy = lo.num_or("dvy", 0);
          l.dvz = lo.num_or("dvz", 0);
          l.dthx = lo.num_or("theta_x", 0);
          if (lo.str_or("frame", "local") == "global") l.frame = LoadFrame::Global;
        } else if (type == "prestrain_linear") {
          l = Load::prestrain_linear(lo.at("bar").as_int(), lo.num_or("strain_i", 0),
                                     lo.num_or("kappa_y_i", 0), lo.num_or("kappa_z_i", 0),
                                     lo.num_or("strain_j", 0), lo.num_or("kappa_y_j", 0),
                                     lo.num_or("kappa_z_j", 0));
          l.kappa_x = lo.num_or("kappa_x_i", 0);
          l.kappa_x2 = lo.num_or("kappa_x_j", 0);
        } else if (type == "prestrain_parabolic") {
          l = Load::prestrain_parabolic(
              lo.at("bar").as_int(),
              Vec3(lo.num_or("strain_i", 0), lo.num_or("kappa_y_i", 0), lo.num_or("kappa_z_i", 0)),
              Vec3(lo.num_or("strain_m", 0), lo.num_or("kappa_y_m", 0), lo.num_or("kappa_z_m", 0)),
              Vec3(lo.num_or("strain_j", 0), lo.num_or("kappa_y_j", 0), lo.num_or("kappa_z_j", 0)));
          l.kappa_x = lo.num_or("kappa_x_i", 0);
          l.kappa_x_m = lo.num_or("kappa_x_m", 0);
          l.kappa_x2 = lo.num_or("kappa_x_j", 0);
        } else if (type == "surface_pressure") {
          std::string fr = lo.str_or("frame", "global");
          SurfaceFrame sf = fr == "local_normal"
                                ? SurfaceFrame::LocalNormal
                                : (fr == "projected" ? SurfaceFrame::Projected : SurfaceFrame::Global);
          l = Load::surface_pressure(lo.at("slab").as_int(), to_vec3(lo.at("w")), sf);
        } else if (type == "shell_prestress") {
          l = Load::shell_prestress(lo.at("tendon").as_int(), lo.num_or("factor", 1.0));
        } else if (type == "shell_imposed") {
          Eigen::Matrix<double, 18, 1> d;
          const auto& a = lo.at("d").as_array();
          require(a.size() == 18, "shell_imposed: 18 values");
          for (int i = 0; i < 18; ++i) d(i) = a[static_cast<std::size_t>(i)].as_number();
          l = Load::shell_imposed(lo.at("shell").as_int(), d);
        } else if (type == "polygon_pressure") {
          std::vector<Vec3> poly;
          for (const auto& p : lo.at("polygon").as_array()) poly.push_back(to_vec3(p));
          std::vector<std::vector<Vec3>> holes;
          if (const json::Value* hs = lo.get("holes"))
            for (const auto& h : hs->as_array()) {
              std::vector<Vec3> r;
              for (const auto& p : h.as_array()) r.push_back(to_vec3(p));
              holes.push_back(r);
            }
          const Vec3 proj = lo.get("proj") ? to_vec3(lo.at("proj")) : Vec3(0, 0, 1);
          l = Load::polygon_pressure(poly, to_vec3(lo.at("q")), holes, lo.get("slab") ? lo.at("slab").as_int() : 0, proj);
        } else if (type == "hydrostatic") {
          l = Load::hydrostatic(lo.num_or("gamma", 0), to_vec3(lo.at("point")), lo.get("slab") ? lo.at("slab").as_int() : 0,
                                lo.num_or("side", 1), lo.get("down") ? to_vec3(lo.at("down")) : Vec3(0, 0, -1));
        } else if (type == "surface_thermal") {
          l = Load::surface_thermal(lo.at("slab").as_int(), lo.num_or("dt", 0),
                                    lo.num_or("grad", 0));
        } else {
          fail("json: unknown load type '" + type + "'");
        }
        l.relative = lo.bool_or("relative", false);
        lc.loads.push_back(l);
      }
      m.cases.push_back(std::move(lc));
    }
  if (const json::Value* cb = v.get("combinations"))
    for (const auto& o : cb->as_array()) {
      Combination c;
      c.id = o.at("id").as_int();
      c.name = o.str_or("name", "");
      for (const auto& t : o.at("terms").as_array())
        c.terms.emplace_back(t.as_array()[0].as_int(), t.as_array()[1].as_number());
      m.combos.push_back(std::move(c));
    }
  auto numbers = [](const json::Value& a, std::size_t n, const char* what) {
    const auto& xs = a.as_array();
    require(xs.size() == n, std::string("json: ") + what + " expects " + std::to_string(n) + " numbers");
    std::vector<double> out;
    for (const auto& x : xs) out.push_back(x.as_number());
    return out;
  };
  auto bools6 = [](const json::Value& a, std::array<bool, 6>& out) {
    const auto& xs = a.as_array();
    require(xs.size() == 6, "json: expected 6 booleans");
    for (std::size_t c = 0; c < 6; ++c) out[c] = xs[c].as_bool();
  };
  if (const json::Value* rl = v.get("rigid_links"))
    for (const auto& o : rl->as_array()) {
      RigidLink r;
      r.id = o.at("id").as_int();
      r.master = o.at("master").as_int();
      for (const auto& x : o.at("slaves").as_array()) r.slaves.push_back(x.as_int());
      if (o.has("dofs")) bools6(o.at("dofs"), r.dofs);
      m.rigid_links.push_back(std::move(r));
    }
  if (const json::Value* ls = v.get("links"))
    for (const auto& o : ls->as_array()) {
      ElasticLink el;
      el.id = o.at("id").as_int();
      el.ni = o.at("ni").as_int();
      el.nj = o.at("nj").as_int();
      const auto k = numbers(o.at("k"), 36, "link k");
      for (int r = 0; r < 6; ++r)
        for (int c = 0; c < 6; ++c) el.k(r, c) = k[static_cast<std::size_t>(6 * r + c)];
      if (o.has("axes")) {
        const auto a = numbers(o.at("axes"), 9, "link axes");
        Mat3 ax;
        for (int r = 0; r < 3; ++r)
          for (int c = 0; c < 3; ++c) ax(r, c) = a[static_cast<std::size_t>(3 * r + c)];
        el.axes = ax;
      }
      el.roll = o.num_or("roll", 0.0);
      el.s = o.num_or("s", 0.5);
      if (const json::Value* laws = o.get("laws")) {
        const auto& ls6 = laws->as_array();
        require(ls6.size() == 6, "json: link laws expects 6 entries");
        for (std::size_t c = 0; c < 6; ++c) {
          if (ls6[c].is_null()) continue;
          LinkLaw law;
          for (const auto& p : ls6[c].as_array())
            law.pts.emplace_back(p.as_array().at(0).as_number(), p.as_array().at(1).as_number());
          el.law[c] = std::move(law);
        }
      }
      m.links.push_back(std::move(el));
    }
  if (const json::Value* cs = v.get("constraints"))
    for (const auto& o : cs->as_array()) {
      Constraint c;
      c.id = o.at("id").as_int();
      for (const auto& t : o.at("terms").as_array()) {
        const auto& x = t.as_array();
        c.terms.push_back(ConstraintTerm{static_cast<Id>(x.at(0).as_int()), static_cast<int>(x.at(1).as_int()),
                                         x.at(2).as_number()});
      }
      m.constraints.push_back(std::move(c));
    }
  if (const json::Value* sk = v.get("skew_supports"))
    for (const auto& o : sk->as_array()) {
      SkewSupport s;
      s.id = o.at("id").as_int();
      s.node = o.at("node").as_int();
      const auto a = numbers(o.at("axes"), 9, "skew support axes");
      for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) s.axes(r, c) = a[static_cast<std::size_t>(3 * r + c)];
      bools6(o.at("fixed"), s.fixed);
      m.skew_supports.push_back(s);
    }
  if (const json::Value* nm = v.get("nodal_masses"))
    for (const auto& o : nm->as_array()) {
      NodalMass n;
      n.node = o.at("node").as_int();
      const auto a = numbers(o.at("m"), 36, "nodal mass");
      for (int r = 0; r < 6; ++r)
        for (int c = 0; c < 6; ++c) n.m(r, c) = a[static_cast<std::size_t>(6 * r + c)];
      n.add = o.bool_or("add", false);
      m.nodal_masses.push_back(n);
    }
  if (const json::Value* bms = v.get("bar_masses"))
    for (const auto& o : bms->as_array()) {
      BarMass b;
      b.bar = o.at("bar").as_int();
      const std::string kind = o.str_or("kind", "auto");
      if (kind == "auto") b.kind = BarMassKind::Auto;
      else if (kind == "nul") b.kind = BarMassKind::Nul;
      else if (kind == "uni") b.kind = BarMassKind::Uni;
      else if (kind == "lin") b.kind = BarMassKind::Lin;
      else if (kind == "parabol") b.kind = BarMassKind::Parabol;
      else fail("json: unknown bar mass kind '" + kind + "'");
      b.m1 = o.num_or("m1", 0.0);
      b.m2 = o.num_or("m2", 0.0);
      b.m3 = o.num_or("m3", 0.0);
      if (o.has("ecc_or")) b.ecc_or = to_vec3(o.at("ecc_or"));
      if (o.has("ecc_ex")) b.ecc_ex = to_vec3(o.at("ecc_ex"));
      const std::string form = o.str_or("form", "model");
      if (form == "model") b.form = BarMass::Form::ModelWide;
      else if (form == "diagonal") b.form = BarMass::Form::Diagonal;
      else if (form == "consistent") b.form = BarMass::Form::Consistent;
      else fail("json: unknown bar mass form '" + form + "'");
      m.bar_masses.push_back(b);
    }
  m.bar_mass_selective = v.get("bar_mass_selective") ? v.at("bar_mass_selective").as_bool() : false;
  return m;
}

const char* const kResultsSignConvention =
    "station forces are the resultant on the face of outward normal +x of the local frame: N > 0 in tension, "
    "My = integral of sigma z dA (My > 0 stretches the +z fibre), dMy/dx = Vz; reactions are the forces of the "
    "supports on the structure, global frame";

namespace {

// ---- result store: one writer and one reader per record ----

template <class Vector>
json::Value vector_json(const Vector& v) {
  json::Value row = json::Value::array();
  for (Eigen::Index c = 0; c < v.size(); ++c) row.push_back(v(c));
  return row;
}

template <class Vector>
Vector vector_from(const json::Value& v, const char* what) {
  const json::Array& row = v.as_array();
  Vector out;
  require(static_cast<Eigen::Index>(row.size()) == out.size(),
          std::string("result store: ") + what + " needs " + std::to_string(out.size()) + " values");
  for (Eigen::Index c = 0; c < out.size(); ++c) out(c) = row[static_cast<std::size_t>(c)].as_number();
  return out;
}

using StationRow = Eigen::Matrix<double, 13, 1>;   // x N Vy Vz T My Mz u v w rx ry rz
using ShellRow = Eigen::Matrix<double, 14, 1>;     // nx ny nxy mx my mxy qx qy m1 m2 m_theta n1 n2 n_theta

json::Value station_json(const Station& s) {
  StationRow row;
  row << s.x, s.N, s.Vy, s.Vz, s.T, s.My, s.Mz, s.defl, s.rot;
  return vector_json(row);
}

Station station_from(const json::Value& v) {
  const StationRow row = vector_from<StationRow>(v, "a station");
  Station s;
  s.x = row(0), s.N = row(1), s.Vy = row(2), s.Vz = row(3), s.T = row(4), s.My = row(5), s.Mz = row(6);
  s.defl = row.segment<3>(7);
  s.rot = row.segment<3>(10);
  return s;
}

json::Value bar_json(const BarResult& b) {
  json::Value o = json::Value::object();
  o.set("bar", b.bar);
  o.set("end_forces", vector_json(b.end_forces_local));
  json::Value stations = json::Value::array(), soil_p = json::Value::array(), soil_m = json::Value::array();
  for (const Station& s : b.stations) stations.push_back(station_json(s));
  for (const Vec3& p : b.soil_p) soil_p.push_back(vector_json(p));
  for (double m : b.soil_m) soil_m.push_back(m);
  o.set("stations", std::move(stations));
  if (!b.soil_p.empty()) o.set("soil_p", std::move(soil_p));
  if (!b.soil_m.empty()) o.set("soil_m", std::move(soil_m));
  return o;
}

BarResult bar_from(const json::Value& v) {
  BarResult b;
  b.bar = static_cast<Id>(v.at("bar").as_int());
  b.end_forces_local = vector_from<Vec12>(v.at("end_forces"), "bar end forces");
  for (const json::Value& s : v.at("stations").as_array()) b.stations.push_back(station_from(s));
  if (const json::Value* p = v.get("soil_p"))
    for (const json::Value& row : p->as_array()) b.soil_p.push_back(vector_from<Vec3>(row, "a soil pressure"));
  if (const json::Value* m = v.get("soil_m"))
    for (const json::Value& x : m->as_array()) b.soil_m.push_back(x.as_number());
  return b;
}

json::Value shell_json(const ShellResult& s) {
  ShellRow row;
  row << s.r.nx, s.r.ny, s.r.nxy, s.r.mx, s.r.my, s.r.mxy, s.r.qx, s.r.qy, s.r.m1, s.r.m2, s.r.m_theta, s.r.n1, s.r.n2,
      s.r.n_theta;
  json::Value o = json::Value::object();
  o.set("shell", s.shell);
  o.set("slab", s.slab);
  o.set("resultants", vector_json(row));
  return o;
}

ShellResult shell_from(const json::Value& v) {
  const ShellRow row = vector_from<ShellRow>(v.at("resultants"), "shell resultants");
  ShellResult s;
  s.shell = static_cast<Id>(v.at("shell").as_int());
  s.slab = static_cast<Id>(v.at("slab").as_int());
  s.r.nx = row(0), s.r.ny = row(1), s.r.nxy = row(2), s.r.mx = row(3), s.r.my = row(4), s.r.mxy = row(5);
  s.r.qx = row(6), s.r.qy = row(7), s.r.m1 = row(8), s.r.m2 = row(9), s.r.m_theta = row(10);
  s.r.n1 = row(11), s.r.n2 = row(12), s.r.n_theta = row(13);
  return s;
}

json::Value case_json(const CaseResult& r) {
  json::Value o = json::Value::object();
  o.set("id", r.id);
  o.set("name", r.name);
  if (r.non_superposable) o.set("non_superposable", true);
  json::Value disp = json::Value::array(), reactions = json::Value::array(), bars = json::Value::array();
  for (const Vec6& d : r.node_disp) disp.push_back(vector_json(d));
  for (const Reaction& re : r.reactions) {
    json::Value ro = json::Value::object();
    ro.set("node", re.node);
    ro.set("R", vector_json(re.R));
    reactions.push_back(std::move(ro));
  }
  for (const BarResult& b : r.bars) bars.push_back(bar_json(b));
  o.set("node_disp", std::move(disp));
  o.set("reactions", std::move(reactions));
  o.set("bars", std::move(bars));
  if (!r.shells.empty()) {
    json::Value shells = json::Value::array();
    for (const ShellResult& s : r.shells) shells.push_back(shell_json(s));
    o.set("shells", std::move(shells));
  }
  if (!r.links.empty()) {
    json::Value links = json::Value::array();
    for (const LinkResult& l : r.links) {
      json::Value lo = json::Value::object();
      lo.set("link", l.link);
      lo.set("deformation", vector_json(l.deformation));
      lo.set("force", vector_json(l.force));
      links.push_back(std::move(lo));
    }
    o.set("links", std::move(links));
  }
  if (!r.constraint_forces.empty()) {
    json::Value forces = json::Value::array();
    for (double f : r.constraint_forces) forces.push_back(f);
    o.set("constraint_forces", std::move(forces));
  }
  if (!r.skew_reactions.empty()) {
    json::Value skew = json::Value::array();
    for (const CaseResult::SkewReaction& s : r.skew_reactions) {
      json::Value so = json::Value::object();
      so.set("node", s.node);
      so.set("local", vector_json(s.local));
      so.set("global", vector_json(s.global));
      skew.push_back(std::move(so));
    }
    o.set("skew_reactions", std::move(skew));
  }
  return o;
}

CaseResult case_from(const json::Value& v) {
  CaseResult r;
  r.id = static_cast<Id>(v.at("id").as_int());
  r.name = v.at("name").as_string();
  r.non_superposable = v.bool_or("non_superposable", false);
  for (const json::Value& d : v.at("node_disp").as_array())
    r.node_disp.push_back(vector_from<Vec6>(d, "a node displacement"));
  for (const json::Value& re : v.at("reactions").as_array())
    r.reactions.push_back({static_cast<Id>(re.at("node").as_int()), vector_from<Vec6>(re.at("R"), "a reaction")});
  for (const json::Value& b : v.at("bars").as_array()) r.bars.push_back(bar_from(b));
  if (const json::Value* shells = v.get("shells"))
    for (const json::Value& s : shells->as_array()) r.shells.push_back(shell_from(s));
  if (const json::Value* links = v.get("links"))
    for (const json::Value& l : links->as_array())
      r.links.push_back({static_cast<Id>(l.at("link").as_int()),
                         vector_from<Vec6>(l.at("deformation"), "a link deformation"),
                         vector_from<Vec6>(l.at("force"), "a link force")});
  if (const json::Value* forces = v.get("constraint_forces"))
    for (const json::Value& f : forces->as_array()) r.constraint_forces.push_back(f.as_number());
  if (const json::Value* skew = v.get("skew_reactions"))
    for (const json::Value& s : skew->as_array())
      r.skew_reactions.push_back({static_cast<Id>(s.at("node").as_int()),
                                  vector_from<Vec6>(s.at("local"), "a skew reaction"),
                                  vector_from<Vec6>(s.at("global"), "a skew reaction")});
  return r;
}

std::string default_generator() {
#ifdef KETA_VERSION
  return std::string("keta ") + KETA_VERSION;
#else
  return "keta";
#endif
}

}  // namespace

json::Value results_to_json(const StaticResults& results, const ResultsProvenance& provenance) {
  json::Value root = json::Value::object();
  root.set("schema_version", kResultsSchemaVersion);
  json::Value origin = json::Value::object();
  origin.set("generator", provenance.generator.empty() ? default_generator() : provenance.generator);
  origin.set("note", provenance.note);
  root.set("provenance", std::move(origin));
  root.set("sign_convention", kResultsSignConvention);
  json::Value cases = json::Value::array(), combos = json::Value::array(), auto_fixed = json::Value::array();
  for (const CaseResult& r : results.cases) cases.push_back(case_json(r));
  for (const CaseResult& r : results.combos) combos.push_back(case_json(r));
  for (const AutoFixedDof& d : results.auto_fixed) {
    json::Value o = json::Value::object();
    o.set("node", d.node);
    o.set("comp", d.comp);
    auto_fixed.push_back(std::move(o));
  }
  root.set("cases", std::move(cases));
  if (!results.combos.empty()) root.set("combinations", std::move(combos));
  if (!results.auto_fixed.empty()) root.set("auto_fixed", std::move(auto_fixed));
  return root;
}

StoredResults results_from_json(const json::Value& v) {
  require(v.is_object() && v.has("schema_version"), "result store: not a keta result document (no schema_version)");
  StoredResults out;
  out.schema_version = static_cast<int>(v.at("schema_version").as_int());
  require(out.schema_version >= 1 && out.schema_version <= kResultsSchemaVersion,
          "result store: schema version " + std::to_string(out.schema_version) + " is not supported (this keta reads up to " +
              std::to_string(kResultsSchemaVersion) + ")");
  out.provenance.generator = v.at("provenance").at("generator").as_string();
  out.provenance.note = v.at("provenance").at("note").as_string();
  out.sign_convention = v.at("sign_convention").as_string();
  for (const json::Value& c : v.at("cases").as_array()) out.results.cases.push_back(case_from(c));
  if (const json::Value* combos = v.get("combinations"))
    for (const json::Value& c : combos->as_array()) out.results.combos.push_back(case_from(c));
  if (const json::Value* fixed = v.get("auto_fixed"))
    for (const json::Value& d : fixed->as_array())
      out.results.auto_fixed.push_back({static_cast<Id>(d.at("node").as_int()), static_cast<int>(d.at("comp").as_int())});
  return out;
}

Model load_model_file(const std::string& path) {
  std::ifstream f(path);
  require(f.good(), "cannot open '" + path + "'");
  std::stringstream ss;
  ss << f.rdbuf();
  return model_from_json(json::Value::parse(ss.str()));
}

StoredResults load_results_file(const std::string& path) {
  std::ifstream f(path);
  require(f.good(), "cannot open '" + path + "'");
  std::stringstream ss;
  ss << f.rdbuf();
  return results_from_json(json::Value::parse(ss.str()));
}

void save_json_file(const std::string& path, const json::Value& v, int indent) {
  std::ofstream f(path);
  require(f.good(), "cannot write '" + path + "'");
  f << v.dump(indent) << "\n";
}

}  // namespace keta
