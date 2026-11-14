// keta Python bindings — variable sections, unilateral and full elastic
// supports, ST1 deck import. Included by keta_py.cpp inside its
// anonymous namespace.
//
//   m.set_tapered(bar, "lin_z", orig=s0, ext=s1)          # CARA VAR LIN Z
//   m.set_support(node, uplift=dict(dz=+1))                # reaction >= 0 only
//   m.set_support(node, decol=("z", +1))                   # all-or-nothing uplift
//   m.set_support(node, stiffness=[kx, ky, kz, 0, 0, 0], rot=math.pi / 4)
//   r = m.pse({bar: dict(zones=[dict(kfz=5e7)], plastic=dict(side="neg", ffz=2e5))})
//   model, info = keta.import_st1(text)

SectionProps props_of(const PyModel& p, const py::handle& h) {
  if (py::isinstance<py::int_>(h)) {
    const Id id = h.cast<Id>();
    for (const Section& s : p.m.sections)
      if (s.id == id) return s.props;
    throw py::key_error("no section " + std::to_string(id));
  }
  const py::dict d = h.cast<py::dict>();
  for (const auto& kv : d) {
    const std::string k = kv.first.cast<std::string>();
    if (k != "A" && k != "Avy" && k != "Avz" && k != "Iy" && k != "Iz" && k != "J")
      throw py::value_error("section properties are A, Avy, Avz, Iy, Iz, J (got '" + k + "')");
  }
  SectionProps s;
  s.A = opt_get<double>(d, "A", 0.0);
  s.Avy = opt_get<double>(d, "Avy", 0.0);
  s.Avz = opt_get<double>(d, "Avz", 0.0);
  s.Iy = opt_get<double>(d, "Iy", 0.0);
  s.Iz = opt_get<double>(d, "Iz", 0.0);
  s.J = opt_get<double>(d, "J", 0.0);
  return s;
}

int dof_index(const std::string& k) {
  static const char* names[6] = {"dx", "dy", "dz", "rx", "ry", "rz"};
  for (int i = 0; i < 6; ++i)
    if (k == names[i]) return i;
  throw py::value_error("support component 'dx', 'dy', 'dz', 'rx', 'ry' or 'rz' (got '" + k + "')");
}

struct PyPse {
  ModelPtr model;
  std::shared_ptr<const PseResults> r;
};
struct PyPseCase {
  ModelPtr model;
  std::shared_ptr<const PseResults> r;
  std::size_t idx = 0;
  const PseCaseResult& get() const { return r->cases[idx]; }
};

SoilZoneMap soil_of(const py::dict& soil) {
  SoilZoneMap out;
  for (const auto& kv : soil) {
    const Id bar = kv.first.cast<Id>();
    const py::dict d = kv.second.cast<py::dict>();
    BarSoil bs;
    for (const py::handle h : d["zones"]) {
      const py::dict z = h.cast<py::dict>();
      SoilZone sz;
      sz.kfx = opt_get<double>(z, "kfx", 0.0);
      sz.kfy = opt_get<double>(z, "kfy", 0.0);
      sz.kfz = opt_get<double>(z, "kfz", 0.0);
      sz.kmx = opt_get<double>(z, "kmx", 0.0);
      sz.xl_end = opt_get<double>(z, "xl", 0.0);
      bs.zones.push_back(sz);
    }
    if (d.contains("plastic") && !d["plastic"].is_none()) {
      const py::dict pl = d["plastic"].cast<py::dict>();
      SoilPlastic sp;
      const std::string side = opt_get<std::string>(pl, "side", "both");
      if (side == "pos") sp.pos = true;
      else if (side == "neg") sp.neg = true;
      else if (side != "both") throw py::value_error("soil side 'pos', 'neg' or 'both'");
      sp.ffx = opt_get<double>(pl, "ffx", 0.0);
      sp.ffy = opt_get<double>(pl, "ffy", 0.0);
      sp.ffz = opt_get<double>(pl, "ffz", 0.0);
      sp.fmx = opt_get<double>(pl, "fmx", 0.0);
      bs.plastic = sp;
    }
    out.emplace(bar, std::move(bs));
  }
  return out;
}

py::dict soil_result(const PseBarSoil& s) {
  py::dict d;
  d["bar"] = s.bar;
  d["bedded"] = s.bedded;
  std::vector<std::array<double, 3>> p;
  for (const Vec3& v : s.soil_p) p.push_back({v.x(), v.y(), v.z()});
  d["p"] = p;
  d["m"] = s.soil_m;
  const Vec3& f = s.total_force_local;
  const Vec3& g = s.total_force_global;
  d["total_local"] = std::array<double, 3>{f.x(), f.y(), f.z()};
  d["total_global"] = std::array<double, 3>{g.x(), g.y(), g.z()};
  d["total_torque"] = s.total_torque_local;
  py::list segs;
  static const char* names[3] = {"elastic", "capped", "inactive"};
  for (const BedSegState& b : s.segments) {
    py::dict q;
    q["x0"] = b.x0;
    q["x1"] = b.x1;
    std::vector<std::string> st;
    for (BedState e : b.state) st.emplace_back(names[static_cast<int>(e)]);
    q["state"] = st;  // (fx, mx, fy, fz)
    q["p_mean"] = b.p_mean;
    q["r_mean"] = b.r_mean;
    segs.append(q);
  }
  d["segments"] = segs;
  return d;
}

py::dict st1_info(const St1Import& r) {
  py::dict d;
  d["title"] = r.title;
  d["warnings"] = r.warnings;
  d["exec_requested"] = r.exec_requested;
  d["exec_modal"] = r.exec_modal;
  d["exec_buckling"] = r.exec_buckling;
  d["etude_pas"] = r.etude_pas;
  return d;
}

void bind_st1(py::module_& mod) {
  auto M = py::reinterpret_borrow<py::class_<PyModel>>(mod.attr("Model"));

  M.def(
      "set_tapered",
      [](PyModel& p, Id bar, const std::string& law, const py::object& orig, const py::object& ext,
         const py::object& mid) {
        TaperedVar v;
        if (law == "lin_y") v.law = TaperLaw::LinY;
        else if (law == "lin_z") v.law = TaperLaw::LinZ;
        else if (law == "para") v.law = TaperLaw::Para;
        else throw py::value_error("tapered law 'lin_y', 'lin_z' or 'para'");
        v.orig = props_of(p, orig);
        v.ext = props_of(p, ext);
        if (v.law == TaperLaw::Para) {
          if (mid.is_none()) throw py::value_error("the 'para' law needs the mid-span properties 'mid'");
          v.mid = props_of(p, mid);
        } else if (!mid.is_none()) {
          throw py::value_error("'mid' is read by the 'para' law only");
        }
        engine([&] { install_tapered(p.m, bar, v); });
      },
      py::arg("bar"), py::arg("law"), py::arg("orig"), py::arg("ext"), py::arg("mid") = py::none(),
      "Variable section along a bar (exact tapered kernel): law 'lin_y' / 'lin_z' (a rectangle whose dimension "
      "along local y / z varies linearly: A, Avy, Avz linear, the inertias cubic or linear accordingly) or "
      "'para' (the inverse of each property is a parabola through orig, mid, ext). orig / mid / ext: a section "
      "id or dict(A, Avy, Avz, Iy, Iz, J). Call once the bar's nodes, roll, eccentricities and material are "
      "final. The bar's own section still feeds its mass and geometric stiffness (section-averaged modal and "
      "buckling analyses); the bar must carry no end release.");

  M.def(
      "set_support",
      [](PyModel& p, Id node, std::optional<py::dict> uplift, std::optional<std::pair<std::string, int>> decol,
         std::optional<std::vector<double>> stiffness, double rot, std::optional<py::dict> cap) {
        Support* s = nullptr;
        for (auto it = p.m.supports.rbegin(); it != p.m.supports.rend(); ++it)
          if (it->node == node) {
            s = &*it;
            break;
          }
        if (!s) {
          Support fresh;
          fresh.node = node;
          p.m.supports.push_back(fresh);
          s = &p.m.supports.back();
        }
        if (uplift) {
          s->oneside = {};
          for (const auto& kv : *uplift) {
            const int c = dof_index(kv.first.cast<std::string>());
            const int sign = kv.second.cast<int>();
            if (sign != 1 && sign != -1 && sign != 0)
              throw py::value_error("uplift sign +1 (reaction >= 0), -1 (reaction <= 0) or 0");
            if (sign != 0 && !s->fixed[static_cast<std::size_t>(c)])
              throw py::value_error("a one-sided component must be blocked by add_support first");
            s->oneside[static_cast<std::size_t>(c)] = sign;
          }
        }
        if (decol) {
          const std::string& ax = decol->first;
          if (ax == "x") s->decol_axis = 0;
          else if (ax == "y") s->decol_axis = 1;
          else if (ax == "z") s->decol_axis = 2;
          else throw py::value_error("decol axis 'x', 'y' or 'z'");
          if (decol->second != 1 && decol->second != -1) throw py::value_error("decol sign +1 or -1");
          s->decol_sign = decol->second;
        }
        if (stiffness) {
          const std::vector<double>& v = *stiffness;
          if (v.size() == 6) s->spring6 = engine([&] { return spring_from_diagonal(v, false); });
          else if (v.size() == 21) s->spring6 = engine([&] { return spring_from_lower_triangle(v, false); });
          else throw py::value_error("stiffness: 6 diagonal terms or the 21 terms of the lower triangle, row by row");
        }
        s->rot = rot;
        if (cap) {
          s->cap = {};
          for (const auto& kv : *cap) {
            const int c = dof_index(kv.first.cast<std::string>());
            const double v = kv.second.cast<double>();
            if (!(v >= 0)) throw py::value_error("a support force cap is >= 0");
            if (v > 0 && !s->fixed[static_cast<std::size_t>(c)])
              throw py::value_error("a force cap needs a component blocked by add_support");
            s->cap[static_cast<std::size_t>(c)] = v;
          }
        }
        p.m.invalidate_index();
      },
      py::arg("node"), py::kw_only(), py::arg("uplift") = py::none(), py::arg("decol") = py::none(),
      py::arg("stiffness") = py::none(), py::arg("rot") = 0.0, py::arg("cap") = py::none(),
      "Completes the support of a node (created if absent). uplift=dict(dz=+1): that blocked component only "
      "pushes (reaction >= 0; -1 for <= 0), solved by active set per load case, non-superposable; "
      "decol=('z', +1): all-or-nothing uplift, every DOF released once the reaction on that global axis changes "
      "sign; stiffness: full elastic support, 6 diagonal terms or 21 lower-triangle terms (dx dy dz rx ry rz) "
      "in the support frame, rotated by `rot` [rad] about global Z; cap=dict(dz=Rmax): the force cap of a blocked "
      "component (bearing yield), past which the support carries exactly +-Rmax while the node moves into it "
      "(on a one-sided component, on its admissible side).");

  M.def(
      "pse",
      [](PyModel& self, const py::dict& soil, int sub_per_zone, int max_iter, double tol, bool exact) {
        ModelPtr model = self.solver_model();
        const SoilZoneMap sm = soil_of(soil);
        PseOptions o;
        o.sub_per_zone = sub_per_zone;
        o.max_iterations = max_iter;
        o.tol = tol;
        o.exact = exact;
        auto r = engine([&] { return std::make_shared<const PseResults>(solve_static_pse(*model, sm, {}, o)); });
        return PyPse{model, r};
      },
      py::arg("soil"), py::arg("sub_per_zone") = 10, py::arg("max_iter") = 50, py::arg("tol") = 1e-8,
      py::arg("exact") = true,
      "Statics of every load case on zoned and pseudo-plastic soils. soil = {bar: dict(zones=[dict(kfx, kfy, "
      "kfz, kmx [Pa, N m/rad per m], xl = end abscissa [m]; the last zone runs to the bar end)], plastic=dict("
      "side 'pos' | 'neg' | 'both', ffx, ffy, ffz [N/m], fmx [N m/m] thresholds))}. The soil reaction "
      "intensity is p = k u_local (sign of the local displacement); side 'pos' keeps the soil where p >= 0, "
      "'neg' where p <= 0. Each zone is cut into sub_per_zone active-set segments. A one-sided or capped soil "
      "is not superposable: no combinations are formed then.");

  py::class_<PyPse>(mod, "PseResult", "Statics on zoned / pseudo-plastic soils, per load case.")
      .def("case",
           [](const PyPse& s, Id id) {
             for (std::size_t i = 0; i < s.r->cases.size(); ++i)
               if (s.r->cases[i].result.id == id) return PyPseCase{s.model, s.r, i};
             throw py::key_error("no load case " + std::to_string(id));
           })
      .def_property_readonly("n_combos", [](const PyPse& s) { return s.r->combos.size(); });
  {
    py::class_<PyPseCase> pc(mod, "PseCase", "One load case on zoned / pseudo-plastic soil.");
    pc.def_property_readonly("iterations", [](const PyPseCase& s) { return s.get().iterations; })
        .def_property_readonly("plastic", [](const PyPseCase& s) { return s.get().plastic; })
        .def("soil",
             [](const PyPseCase& s, Id bar) {
               for (const PseBarSoil& b : s.get().soil)
                 if (b.bar == bar) return soil_result(b);
               throw py::key_error("no bar " + std::to_string(bar));
             },
             "dict(p per station (x, y, z local) [N/m], m [N m/m], total_local, total_global (integral of p; the "
             "soil pushes the structure by its negative), total_torque, segments=[dict(x0, x1, state of (fx, mx, "
             "fy, fz), p_mean, r_mean)]).");
    bind_case_accessors(
        pc, [](const PyPseCase& s) { return std::pair<ModelPtr, const CaseResult&>(s.model, s.get().result); },
        "bar_end_forces");
  }

  mod.def(
      "import_st1",
      [](const std::string& text, double gravity, bool is_path) {
        St1Options o;
        o.gravity = gravity;
        St1Import r = engine([&] { return is_path ? import_st1_file(text, o) : import_st1(text, o); });
        return py::make_tuple(PyModel::from_model(r.model), st1_info(r));
      },
      py::arg("text"), py::arg("gravity") = 1.0, py::arg("is_path") = false,
      "Runs an ST1 deck (text, or a file path with is_path=True; LIRE resolved next to it) and returns "
      "(Model, info): the model with every CHARG as a load case and every COMB as a combination, info = "
      "dict(title, warnings, exec_requested, exec_modal, exec_buckling, etude_pas). gravity scales POIDS PROPRE.");

  mod.def(
      "run_st1",
      [](const std::string& text, bool st1_edition, double gravity) {
        St1Options o;
        o.gravity = gravity;
        o.st1_edition = st1_edition;
        return engine([&] {
          st1::St1Session s(o);
          s.run(text);
          return s.runtime->output().console;
        });
      },
      py::arg("text"), py::arg("st1_edition") = false, py::arg("gravity") = 1.0,
      "Runs an ST1 deck and returns its console listing (RESU, GET, LISTER ...); st1_edition prints the tables "
      "in ST1's own layout and sign conventions.");
}
