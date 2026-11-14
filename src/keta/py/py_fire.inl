// keta Python bindings — fire: EN 1991-1-2 curves, EN 1993-1-2 /
// EN 1992-1-2 reduction factors and thermal strains, steel member heating,
// and the frame fire analysis (section heat conduction, hot rigidities and
// thermal strains, frame iteration). Included by keta_py.cpp inside its
// anonymous namespace.
//
// Section frame: the thickness runs along the bar's local y (face 'sup' = +y,
// 'inf' = -y), the width along local z (faces 'pz', 'mz'). A beam bending in
// a vertical plane is rolled so its local y is vertical: beta = pi/2 on a
// horizontal bar puts local y on +Z, 'sup' on top.
//
//   b = m.add_beam(n1, n2, section=s, material=c, beta=math.pi / 2)
//   r = m.fire(bars=[dict(bar=b, section=dict(kind="rectangle", thickness=0.6,
//               width=0.3, exposures=dict(inf=0, pz=0, mz=0)), fck=30e6,
//               rebars=[dict(diameter=0.02, count=3, placement="inf", cover=0.04)])],
//              exposures=[dict(kind="iso834")], duration=90, preload_case=1)

struct PyFire {
  ModelPtr model;
  std::shared_ptr<const FireOutcome> r;
};

ExposureCurve curve_of(const std::string& s, const std::vector<std::pair<double, double>>& table) {
  ExposureCurve c;
  if (s == "iso834") c.kind = ExposureCurve::Kind::Iso834;
  else if (s == "external") c.kind = ExposureCurve::Kind::External;
  else if (s == "hydrocarbon") c.kind = ExposureCurve::Kind::Hydrocarbon;
  else if (s == "hcm") c.kind = ExposureCurve::Kind::Hcm;
  else if (s == "table") {
    c.kind = ExposureCurve::Kind::Table;
    c.table = table;
  } else
    throw py::value_error("fire curve 'iso834', 'external', 'hydrocarbon', 'hcm' or 'table'");
  return c;
}
Aggregate aggregate_of(const std::string& s) {
  if (s == "siliceous") return Aggregate::Siliceous;
  if (s == "calcareous") return Aggregate::Calcareous;
  throw py::value_error("aggregate 'siliceous' or 'calcareous'");
}
SteelClass steel_class_of(const std::string& s) {
  if (s == "N") return SteelClass::N;
  if (s == "X") return SteelClass::X;
  throw py::value_error("reinforcement class 'N' or 'X'");
}
SteelForming forming_of(const std::string& s) {
  if (s == "hot") return SteelForming::HotRolled;
  if (s == "cold") return SteelForming::ColdWorked;
  throw py::value_error("forming 'hot' or 'cold'");
}

Exposure exposure_of(const py::dict& d) {
  Exposure e;
  const std::string k = opt_get<std::string>(d, "kind", "iso834");
  if (k == "iso834") e.kind = Exposure::Kind::Standard;
  else if (k == "hydrocarbon") e.kind = Exposure::Kind::Hydrocarbon;
  else if (k == "hcm") e.kind = Exposure::Kind::Hcm;
  else if (k == "adiabatic") e.kind = Exposure::Kind::Adiabatic;
  else if (k == "constant") e.kind = Exposure::Kind::ConstantGas;
  else if (k == "gas_table") e.kind = Exposure::Kind::GasTable;
  else if (k == "wall_table") e.kind = Exposure::Kind::WallTable;
  else throw py::value_error("exposure kind 'iso834', 'hydrocarbon', 'hcm', 'adiabatic', 'constant', 'gas_table' or 'wall_table'");
  e.temp = opt_get<double>(d, "temp", 20.0);
  e.table = opt_get<std::vector<std::pair<double, double>>>(d, "table", {});
  e.protection = opt_get<double>(d, "protection", 0.0);
  if (d.contains("convection")) {
    e.user_convection = true;
    e.alpha_c_user = d["convection"].cast<double>();
  }
  return e;
}

FireBarSpec fire_bar_of(const py::dict& d) {
  FireBarSpec b;
  b.bar = d["bar"].cast<Id>();
  const py::dict s = d["section"].cast<py::dict>();
  const std::string kind = opt_get<std::string>(s, "kind", "rectangle");
  FireSection& f = b.section;
  if (kind == "slab") f.kind = FireSection::Kind::Slab;
  else if (kind == "rectangle") f.kind = FireSection::Kind::Rectangle;
  else if (kind == "perimeter") f.kind = FireSection::Kind::Perimeter;
  else throw py::value_error("fire section kind 'slab', 'rectangle' or 'perimeter'");
  f.thickness = opt_get<double>(s, "thickness", 0.0);
  f.width = opt_get<double>(s, "width", 1.0);
  f.maille = opt_get<double>(s, "mesh", 0.01);
  if (s.contains("exposures")) {
    const py::dict e = s["exposures"].cast<py::dict>();
    f.expo_vy = opt_get<int>(e, "sup", -1);
    f.expo_wy = opt_get<int>(e, "inf", -1);
    for (const char* k : {"py", "my"})
      if (e.contains(k))
        throw py::value_error("fire section faces are 'sup', 'inf' (+y / -y, the thickness) and 'pz', 'mz' "
                              "(+z / -z, the width)");
    f.expo_vz = opt_get<int>(e, "pz", -1);
    f.expo_wz = opt_get<int>(e, "mz", -1);
  }
  if (s.contains("polygon"))
    for (const auto& [y, z] : s["polygon"].cast<std::vector<std::pair<double, double>>>()) f.polygon.push_back(Vec3(y, z, 0));
  f.seg_exposure = opt_get<std::vector<int>>(s, "segment_exposures", {});
  b.material.fck = opt_get<double>(d, "fck", 30e6);
  b.material.tension = opt_get<bool>(d, "tension", false);
  b.material.thermal.aggregate = aggregate_of(opt_get<std::string>(d, "aggregate", "siliceous"));
  b.material.thermal.uu = opt_get<double>(d, "moisture", 0.0);
  b.material.thermal.lambda_alpha = opt_get<double>(d, "lambda_alpha", 1.0);
  if (d.contains("rebars"))
    for (const py::handle h : d["rebars"]) {
      const py::dict r = h.cast<py::dict>();
      RebarLayer L;
      L.steel.diameter = opt_get<double>(r, "diameter", 0.0);
      L.steel.area = opt_get<double>(r, "area", 0.0);
      L.steel.E = opt_get<double>(r, "E", 200e9);
      L.steel.fyk = opt_get<double>(r, "fyk", 500e6);
      L.steel.cls = steel_class_of(opt_get<std::string>(r, "cls", "N"));
      L.steel.forming = forming_of(opt_get<std::string>(r, "forming", "cold"));
      L.count = opt_get<int>(r, "count", 1);
      const std::string pl = opt_get<std::string>(r, "placement", "inf");
      if (pl == "sup") L.placement = RebarLayer::Placement::CoverSup;
      else if (pl == "inf") L.placement = RebarLayer::Placement::CoverInf;
      else if (pl == "local") L.placement = RebarLayer::Placement::Local;
      else throw py::value_error("rebar placement 'sup', 'inf' or 'local'");
      L.cover = opt_get<double>(r, "cover", 0.0);
      L.y = opt_get<double>(r, "y", 0.0);
      L.z = opt_get<double>(r, "z", 0.0);
      b.rebars.push_back(L);
    }
  return b;
}

void bind_fire(py::module_& mod) {
  auto cls = py::reinterpret_borrow<py::class_<PyModel>>(mod.attr("Model"));
  mod.def("fire_gas_temperature",
          [](const std::string& curve, double t, const std::vector<std::pair<double, double>>& table) {
            return curve_of(curve, table).gas_temp(t);
          },
          py::arg("curve"), py::arg("t_minutes"), py::arg("table") = std::vector<std::pair<double, double>>{},
          "Gas temperature [degC] of an EN 1991-1-2 curve at t minutes.");
  mod.def("steel_temperature_unprotected",
          [](const std::string& curve, double am_v, double ksh, double t_end, double dt) {
            return engine([&] { return steel_temperature_unprotected(curve_of(curve, {}), am_v, ksh, t_end, dt); });
          },
          py::arg("curve"), py::arg("am_v"), py::arg("ksh") = 1.0, py::arg("t_end_minutes") = 60.0, py::arg("dt") = 5.0,
          "EN 1993-1-2 §4.2.5.1: steel temperature history every dt seconds.");
  mod.def("steel_temperature_protected",
          [](const std::string& curve, double ap_v, double dp, double lp, double rp, double cp, double t_end, double dt) {
            return engine([&] { return steel_temperature_protected(curve_of(curve, {}), ap_v, dp, lp, rp, cp, t_end, dt); });
          },
          py::arg("curve"), py::arg("ap_v"), py::arg("dp"), py::arg("lambda_p"), py::arg("rho_p"), py::arg("cp_p"),
          py::arg("t_end_minutes") = 60.0, py::arg("dt") = 5.0, "EN 1993-1-2 §4.2.5.2 protected steel.");
  mod.def("steel_ky", &steel_ky, py::arg("theta"));
  mod.def("steel_kp", &steel_kp, py::arg("theta"));
  mod.def("steel_ke", &steel_ke, py::arg("theta"));
  mod.def("steel_critical_temperature", &steel_critical_temperature, py::arg("mu0"));
  mod.def("concrete_kc", [](double t, const std::string& a) { return concrete_kc(t, aggregate_of(a)); }, py::arg("theta"),
          py::arg("aggregate") = "siliceous");
  mod.def("concrete_eps_c1", &concrete_eps_c1, py::arg("theta"));
  mod.def("concrete_eps_cu1", &concrete_eps_cu1, py::arg("theta"));
  mod.def("rebar_ks",
          [](double t, const std::string& c, const std::string& f) { return rebar_ks(t, steel_class_of(c), forming_of(f)); },
          py::arg("theta"), py::arg("cls") = "N", py::arg("forming") = "cold");
  mod.def("rebar_kE",
          [](double t, const std::string& c, const std::string& f) { return rebar_kE(t, steel_class_of(c), forming_of(f)); },
          py::arg("theta"), py::arg("cls") = "N", py::arg("forming") = "cold");
  mod.def("concrete_thermal_strain", [](double t, const std::string& a) { return concrete_thermal_strain(t, aggregate_of(a)); },
          py::arg("theta"), py::arg("aggregate") = "siliceous");
  mod.def("rebar_thermal_strain", &rebar_thermal_strain, py::arg("theta"));

  cls.def(
      "fire",
      [](const PyModel& self, const py::list& bars, const py::list& exposures, double duration, double t_init,
         double dt, std::optional<Id> preload, bool local, double tol, int max_iter, bool centroid_var) {
        FireProblem p;
        for (const py::handle h : exposures) p.exposures.push_back(exposure_of(h.cast<py::dict>()));
        for (const py::handle h : bars) p.bars.push_back(fire_bar_of(h.cast<py::dict>()));
        p.duration = duration;
        p.t_init = t_init;
        p.dt = dt;
        p.preload_case = preload.value_or(0);
        p.local = local;
        FireAnalysisOptions o;
        o.tol = tol;
        o.max_iter = max_iter;
        o.centroid_var = centroid_var;
        ModelPtr model = self.solver_model();
        auto r = engine([&] {
          py::gil_scoped_release nogil;
          return std::make_shared<const FireOutcome>(solve_fire(*model, p, o));
        });
        return PyFire{model, r};
      },
      py::arg("bars"), py::arg("exposures"), py::arg("duration") = 60.0, py::arg("t_init") = 20.0, py::arg("dt") = 30.0,
      py::arg("preload_case") = py::none(), py::arg("local") = false, py::arg("tol") = 1e-3, py::arg("max_iter") = 500,
      py::arg("centroid_var") = true,
      "Fire analysis of the frame: heat conduction in each fire bar's section (EN 1992-1-2 properties), hot "
      "secant rigidities and thermal strains, frame re-solved under the preload case until the displacements "
      "settle (non-convergence = failure). bars: [dict(bar, section=dict(kind 'slab'|'rectangle'|'perimeter', "
      "thickness, width, exposures=dict(sup, inf, pz, mz), polygon, segment_exposures, mesh), fck, tension, "
      "aggregate, moisture, lambda_alpha, rebars=[dict(diameter | area, count, placement 'sup'|'inf'|'local', "
      "cover, y, z, fyk, E, cls, forming)])]; exposures: [dict(kind 'iso834'|'hydrocarbon'|'hcm'|'adiabatic'|"
      "'constant'|'gas_table'|'wall_table', temp, table, protection, convection)]. duration [min]. Section frame: "
      "the thickness runs along the bar's local y (sup = +y), the width along local z (pz = +z); roll a "
      "horizontal beam by beta = pi/2 to put its thickness vertical.");

  py::class_<PyFire>(mod, "FireResult", "Hot equilibrium of a fire analysis.")
      .def("displacement", [](const PyFire& s, Id node) { return dof6(s.r->result.node_disp.at(node_pos(*s.model, node))); })
      .def("bar_forces", [](const PyFire& s, Id bar) { return stations_of(bar_of(s.r->result, bar)); })
      .def("reactions", [](const PyFire& s) { return reactions_of(s.r->result); })
      .def_property_readonly("iterations", [](const PyFire& s) { return s.r->iterations; })
      .def("sections",
           [](const PyFire& s) {
             py::list out;
             for (const HotSectionState& h : s.r->states) {
               py::dict d;
               d["bar"] = h.bar;
               d["EA"] = h.EA;
               d["EIy"] = h.EIy;
               d["EIz"] = h.EIz;
               d["yc"] = h.yc;
               d["zc"] = h.zc;
               d["N_th"] = h.N_th;
               d["M_th_y"] = h.M_th_y;
               d["M_th_z"] = h.M_th_z;
               d["theta_max"] = h.theta_max;
               d["theta_mean"] = h.theta_mean;
               d["theta_rebar_max"] = h.theta_rebar_max;
               out.append(d);
             }
             return out;
           },
           "Hot section state per fire bar: rigidities, stiffness centroid, thermal resultants, temperatures.");
}
