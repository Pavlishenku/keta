// keta Python bindings — staged construction, prestressing tendons and the
// time-dependent concrete laws. Included by keta_py.cpp inside its
// anonymous namespace: it reuses PyModel and the result helpers there.
//
//   m.set_concrete_law(mat, "EC2", fc=35e6, rh=70, h_mm=300, cement="N")
//   t = m.add_tendon([b1, b2, b3], P0=3e6, profile=[(0, 0, -0.4), (30, 0, 0.3)],
//                    Ap=2.8e-3, fpk=1860e6)
//   r = m.phased([
//       dict(name="coulage", time=28, cases=[g]),
//       dict(name="mise en tension", time=35, tension=[dict(tendon=t)]),
//       dict(name="10 000 j", time=10000)])
//   r.phase(2).displacement(node), r.phase(2).tendons()

struct PyPhased {
  ModelPtr model;
  std::shared_ptr<const PhasedResults> r;
};
struct PyPhase {
  ModelPtr model;
  std::shared_ptr<const PhasedResults> r;
  std::size_t idx = 0;
  const PhaseResult& get() const { return r->phases[idx]; }
};

MatLaw::Code law_code(const std::string& s) {
  static const std::map<std::string, MatLaw::Code> codes = {
      {"DIRECT", MatLaw::Code::Direct},     {"BPEL", MatLaw::Code::BPEL},         {"BHP1", MatLaw::Code::BHP1},
      {"BHP2", MatLaw::Code::BHP2},         {"BPEL_BL", MatLaw::Code::BPEL_BL},   {"CEB", MatLaw::Code::CEB},
      {"EC2", MatLaw::Code::EC2},           {"EC2_BHP1", MatLaw::Code::EC2_BHP1}, {"EC2_BHP2", MatLaw::Code::EC2_BHP2},
      {"FIB2010", MatLaw::Code::FIB2010},   {"ACI_B3", MatLaw::Code::ACI_B3},     {"EC2_2027", MatLaw::Code::EC2_2027},
      {"IP2", MatLaw::Code::IP2},           {"MC78", MatLaw::Code::MC78},         {"COURBON", MatLaw::Code::COURBON}};
  std::string u = s;
  for (char& c : u) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  const auto it = codes.find(u);
  if (it == codes.end()) throw py::value_error("unknown concrete law '" + s + "'");
  return it->second;
}

MatLaw make_law(const std::string& code, const py::kwargs& kw) {
  MatLaw L;
  L.code = law_code(code);
  for (const auto& [k, v] : kw) {
    const std::string key = py::str(k);
    if (key == "fc") L.fc = v.cast<double>();
    else if (key == "e_direct") L.e_direct = v.cast<double>();
    else if (key == "rh") L.rh = v.cast<double>();
    else if (key == "rm_cm") L.rm_cm = v.cast<double>();
    else if (key == "h_mm") L.h_mm = v.cast<double>();
    else if (key == "ros") L.ros = v.cast<double>();
    else if (key == "ts_dry") L.ts_dry = v.cast<double>();
    else if (key == "es") L.es = v.cast<double>();
    else if (key == "er") L.er = v.cast<double>();
    else if (key == "theta") L.theta = v.cast<double>();
    else if (key == "kb") L.kb = v.cast<double>();
    else if (key == "k1_consistency") L.k1_consistency = v.cast<double>();
    else if (key == "phi_f2") L.phi_f2 = v.cast<double>();
    else if (key == "silica") L.silica = v.cast<bool>();
    else if (key == "fib_secant") L.fib_secant = v.cast<bool>();
    // affinity coefficients of the law's terms (1 = the law as published)
    else if (key == "k_modulus") L.avance.k_module = v.cast<double>();
    else if (key == "k_shrinkage_autogenous") L.avance.k_retrait_endogene = v.cast<double>();
    else if (key == "k_shrinkage_drying") L.avance.k_retrait_dessication = v.cast<double>();
    else if (key == "k_creep_basic") L.avance.k_fluage_endogene = v.cast<double>();
    else if (key == "k_creep_drying") L.avance.k_fluage_dessication = v.cast<double>();
    else if (key == "k_creep_kinetics") L.avance.k_fluage_cinetique = v.cast<double>();
    else if (key == "k_shrinkage_autogenous_kinetics") L.avance.k_retrait_endogene_cin = v.cast<double>();
    else if (key == "k_shrinkage_drying_kinetics") L.avance.k_retrait_dessication_cin = v.cast<double>();
    else if (key == "k_creep_basic_kinetics") L.avance.k_fluage_endogene_cin = v.cast<double>();
    else if (key == "k_creep_drying_kinetics") L.avance.k_fluage_dessication_cin = v.cast<double>();
    else if (key == "cement") {
      const std::string c = py::str(v);
      if (c == "S") L.cement = CementClass::S;
      else if (c == "N") L.cement = CementClass::N;
      else if (c == "R") L.cement = CementClass::R;
      else if (c == "RS") L.cement = CementClass::RS;
      else throw py::value_error("cement class S, N, R or RS");
    } else if (key != "creep_tension" && key != "shrinkage_law" && key != "ts" && key != "nonlinear")
      throw py::value_error("unknown concrete-law parameter '" + key + "'");
  }
  return L;
}

// a piecewise-linear eccentricity profile (s, e_y, e_z) as polynomial segments
std::vector<PolySegment> linear_profile(const std::vector<std::array<double, 3>>& pts) {
  if (pts.size() < 1) throw py::value_error("tendon profile: at least one point (s, e_y, e_z)");
  std::vector<PolySegment> out;
  if (pts.size() == 1) {
    PolySegment p;
    p.s0 = 0;
    p.s1 = 1e30;
    p.cy[0] = pts[0][1];
    p.cz[0] = pts[0][2];
    return {p};
  }
  for (std::size_t i = 0; i + 1 < pts.size(); ++i) {
    const double ds = pts[i + 1][0] - pts[i][0];
    if (!(ds > 0)) throw py::value_error("tendon profile: s strictly increasing");
    PolySegment p;
    p.s0 = pts[i][0];
    p.s1 = pts[i + 1][0];
    p.cy = {pts[i][1], (pts[i + 1][1] - pts[i][1]) / ds, 0, 0};
    p.cz = {pts[i][2], (pts[i + 1][2] - pts[i][2]) / ds, 0, 0};
    out.push_back(p);
  }
  return out;
}

template <class T>
T opt_get(const py::dict& d, const char* k, T def) {
  return d.contains(k) ? d[k].cast<T>() : def;
}

Phase make_phase(const PyModel& self, const py::dict& d, std::size_t index) {
  static const std::set<std::string> known = {
      "name", "time", "activate", "deactivate", "add_supports", "remove_supports", "cases", "tension", "inject",
      "inject_partial", "activate_age", "regime", "imperfection", "cable_actions", "modify_sections", "place_hinges",
      "close_hinges"};
  for (const auto& [k, v] : d)
    if (!known.count(py::str(k))) throw py::value_error("phase " + std::to_string(index) + ": unknown key '" + std::string(py::str(k)) + "'");
  Phase p;
  p.name = opt_get<std::string>(d, "name", "phase " + std::to_string(index + 1));
  p.time = opt_get<double>(d, "time", 0.0);
  p.activate_bars = opt_get<std::vector<Id>>(d, "activate", {});
  p.deactivate_bars = opt_get<std::vector<Id>>(d, "deactivate", {});
  p.remove_support_nodes = opt_get<std::vector<Id>>(d, "remove_supports", {});
  p.cases = opt_get<std::vector<Id>>(d, "cases", {});
  if (d.contains("add_supports"))
    for (const py::handle h : d["add_supports"]) {
      const py::dict s = h.cast<py::dict>();
      Support sp;
      sp.node = s["node"].cast<Id>();
      const char* names[6] = {"dx", "dy", "dz", "rx", "ry", "rz"};
      for (int c = 0; c < 6; ++c) sp.fixed[static_cast<std::size_t>(c)] = opt_get<bool>(s, names[c], c < 3);
      p.add_supports.push_back(sp);
    }
  if (d.contains("tension"))
    for (const py::handle h : d["tension"]) {
      const py::dict t = h.cast<py::dict>();
      TensionAction a;
      a.tendon = t["tendon"].cast<Id>();
      a.fraction = opt_get<double>(t, "fraction", 1.0);
      const std::string mode = opt_get<std::string>(t, "mode", "jack");
      if (mode == "jack") a.mode = TensionAction::Mode::Jack;
      else if (mode == "bed") a.mode = TensionAction::Mode::Bed;
      else if (mode == "release") a.mode = TensionAction::Mode::Release;
      else throw py::value_error("tension mode 'jack', 'bed' or 'release'");
      a.unbonded = opt_get<bool>(t, "unbonded", false);
      a.group = opt_get<int>(t, "group", 0);
      p.tension.push_back(a);
    }
  p.inject = opt_get<std::vector<Id>>(d, "inject", {});
  if (d.contains("inject_partial"))
    for (const py::handle h : d["inject_partial"]) {
      const auto [t, s0, s1] = h.cast<std::tuple<Id, double, double>>();
      p.inject_partial.push_back({t, s0, s1});
    }
  if (d.contains("activate_age"))
    for (const auto& [k, v] : d["activate_age"].cast<py::dict>()) p.activate_age[k.cast<Id>()] = v.cast<double>();
  const std::string regime = opt_get<std::string>(d, "regime", "linear");
  if (regime == "linear") p.regime = Phase::Regime::Linear;
  else if (regime == "second_order") p.regime = Phase::Regime::SecondOrder;
  else throw py::value_error("regime 'linear' or 'second_order'");
  if (d.contains("imperfection"))
    for (const py::handle h : d["imperfection"]) {
      const auto [n, v] = h.cast<std::pair<Id, std::array<double, 6>>>();
      Vec6 x;
      for (int c = 0; c < 6; ++c) x(c) = v[static_cast<std::size_t>(c)];
      p.imperfection.push_back({n, x});
    }
  if (d.contains("modify_sections"))
    for (const py::handle h : d["modify_sections"]) p.modify_sections.push_back(h.cast<std::pair<Id, Id>>());
  auto hinges = [&](const char* key, std::vector<Phase::HingeAction>& out) {
    if (!d.contains(key)) return;
    for (const py::handle h : d[key]) {
      const py::dict x = h.cast<py::dict>();
      Phase::HingeAction a;
      a.bar = x["bar"].cast<Id>();
      a.end = opt_get<int>(x, "end", 1);
      const std::string dofs = opt_get<std::string>(x, "dofs", "rz");
      for (const char* name : {"dx", "dy", "dz", "rx", "ry", "rz"}) {
        const int c = static_cast<int>(std::string("dxdydzrxryrz").find(name)) / 2;
        if (dofs.find(name) != std::string::npos) a.dofs[static_cast<std::size_t>(c)] = true;
      }
      out.push_back(a);
    }
  };
  hinges("place_hinges", p.place_hinges);
  hinges("close_hinges", p.close_hinges);
  if (d.contains("cable_actions"))
    for (const py::handle h : d["cable_actions"]) {
      const py::dict x = h.cast<py::dict>();
      Phase::CableAction a;
      const Id cable = x["cable"].cast<Id>();
      a.bar = self.cable_bar_id(cable);
      const std::string kind = x["kind"].cast<std::string>();
      using K = Phase::CableAction::Kind;
      if (kind == "length") a.kind = K::Length;
      else if (kind == "horizontal") a.kind = K::Horizontal;
      else if (kind == "end_tension") a.kind = K::EndTension;
      else if (kind == "thermal") a.kind = K::Thermal;
      else if (kind == "rupture") a.kind = K::Rupture;
      else throw py::value_error("cable action: length, horizontal, end_tension, thermal or rupture");
      a.value = opt_get<double>(x, "value", 0.0);
      p.cable_actions.push_back(a);
    }
  return p;
}

void bind_phasing(py::module_& mod) {
  auto cls = py::reinterpret_borrow<py::class_<PyModel>>(mod.attr("Model"));
  cls.def(
         "set_concrete_law",
         [](PyModel& self, Id material, const std::string& code, const py::kwargs& kw) {
           const MatLaw L = make_law(code, kw);
           MatLaw shr = L;
           // the shrinkage law takes the same parameters (fc, rh, h_mm, cement ...)
           if (kw.contains("shrinkage_law")) shr = make_law(kw["shrinkage_law"].cast<std::string>(), kw);
           const double ts = kw.contains("ts") ? kw["ts"].cast<double>() : 0.0;
           const bool nonlinear = kw.contains("nonlinear") && kw["nonlinear"].cast<bool>();
           std::shared_ptr<ConcreteLaw> cl = engine([&] { return concrete_law_from(L, shr, ts, 1.0, nonlinear); });
           if (kw.contains("creep_tension")) cl->creep_tension = kw["creep_tension"].cast<bool>();
           for (Material& x : self.m.materials)
             if (x.id == material) {
               x.concrete_law = cl;
               return;
             }
           throw py::key_error("no material " + std::to_string(material));
         },
         py::arg("material"), py::arg("code") = "EC2",
         "Give a material a time-dependent concrete law (creep phi(t, t0), modulus, shrinkage) for the phased "
         "analysis. code: EC2, EC2_2027, EC2_BHP1/2, CEB, FIB2010, BPEL, BHP1/2, BPEL_BL, ACI_B3, IP2, MC78, "
         "COURBON, DIRECT. Keywords (SI, ages in days): fc, rh [%], h_mm (2Ac/u), rm_cm (Ac/u, BPEL family), "
         "ts_dry, cement ('S','N','R','RS'), e_direct, er, ros, theta, kb, k1_consistency, phi_f2 (MC78), "
         "silica, creep_tension (False: no creep in tension), shrinkage_law (another code for the shrinkage, same "
         "parameters), ts, nonlinear (True: EN 1992-1-1 §3.1.4(4) non-linear creep above 0.45 fcm(t0); default "
         "linear). Affinity coefficients (default 1, multiply the law's terms): k_modulus, k_shrinkage_autogenous, "
         "k_shrinkage_drying, k_creep_basic, k_creep_drying, k_creep_kinetics, and the kinetics factors "
         "k_shrinkage_autogenous_kinetics, k_shrinkage_drying_kinetics, k_creep_basic_kinetics, "
         "k_creep_drying_kinetics (they divide the age in their term's time function, EN 1992-1-1:2027 B.4.10).")
      .def(
          "set_concrete_law_functions",
          [](PyModel& self, Id material, py::function phi, py::function e_ref, std::optional<py::function> shrinkage) {
            auto cl = std::make_shared<ConcreteLaw>();
            cl->phi = [phi](double t, double t0) {
              py::gil_scoped_acquire g;
              return phi(t, t0).cast<double>();
            };
            cl->E_ref = [e_ref](double t0) {
              py::gil_scoped_acquire g;
              return e_ref(t0).cast<double>();
            };
            if (shrinkage) {
              py::function sh = *shrinkage;
              cl->shrinkage = [sh](double t) {
                py::gil_scoped_acquire g;
                return sh(t).cast<double>();
              };
            } else {
              cl->shrinkage = [](double) { return 0.0; };
            }
            for (Material& x : self.m.materials)
              if (x.id == material) {
                x.concrete_law = cl;
                return;
              }
            throw py::key_error("no material " + std::to_string(material));
          },
          py::arg("material"), py::arg("phi"), py::arg("E_ref"), py::arg("shrinkage") = py::none(),
          "A user concrete law as Python functions: phi(t, t0), E_ref(t0) [Pa], shrinkage(t) (< 0 = shortening).")
      .def(
          "add_tendon",
          [](PyModel& self, std::vector<Id> bars, double P0, std::vector<std::array<double, 3>> profile,
             const std::string& jack, double mu, double k, double slip, double Ep, double Ap, double fpk,
             double rho1000, int relaxation_class, double s_begin, double s_end, const std::string& name,
             std::optional<py::list> segments, bool exact_angle, double phi_abs, int non_simul,
             const std::string& relaxation, double relax_teb_days, double thermal_loss,
             std::vector<std::pair<double, double>> free_runs, std::vector<std::pair<double, double>> force_profile) {
            Tendon t;
            t.id = 1;
            for (const Tendon& x : self.m.tendons) t.id = std::max(t.id, x.id + 1);
            t.name = name;
            t.bars = std::move(bars);
            if (t.bars.empty()) throw py::value_error("tendon: at least one bar");
            for (Id b : t.bars) (void)self.m.bar(b);
            if (segments) {
              if (!profile.empty()) throw py::value_error("tendon: give profile= or segments=, not both");
              for (const py::handle& h : *segments) {
                const py::dict d = py::reinterpret_borrow<py::dict>(h);
                PolySegment ps;
                ps.s0 = d["s0"].cast<double>();
                ps.s1 = d["s1"].cast<double>();
                if (!(ps.s1 > ps.s0)) throw py::value_error("tendon segment: s1 > s0");
                if (d.contains("cy")) {
                  const auto c = d["cy"].cast<std::vector<double>>();
                  if (c.size() > 4) throw py::value_error("tendon segment: at most 4 coefficients (cubic)");
                  for (std::size_t i = 0; i < c.size(); ++i) ps.cy[i] = c[i];
                }
                if (d.contains("cz")) {
                  const auto c = d["cz"].cast<std::vector<double>>();
                  if (c.size() > 4) throw py::value_error("tendon segment: at most 4 coefficients (cubic)");
                  for (std::size_t i = 0; i < c.size(); ++i) ps.cz[i] = c[i];
                }
                if (!t.profile.empty() && std::abs(t.profile.back().s1 - ps.s0) > 1e-9)
                  throw py::value_error("tendon segments: contiguous (s0 = previous s1)");
                t.profile.push_back(ps);
              }
              if (t.profile.empty()) throw py::value_error("tendon: segments= is empty");
            } else {
              t.profile = linear_profile(profile);
            }
            t.exact_angle = exact_angle;
            t.losses.phi_abs = phi_abs;
            if (non_simul < 1) throw py::value_error("tendon: non_simul >= 1");
            t.losses.non_simul = non_simul;
            if (relaxation == "bpel") t.losses.relaxation_bpel = true;
            else if (relaxation != "ec2") throw py::value_error("tendon relaxation 'ec2' or 'bpel'");
            t.losses.relax_teb_days = relax_teb_days;
            t.losses.thermal_loss = thermal_loss;
            t.jack_force = P0;
            if (jack == "both") t.jack_both_ends = true;
            else if (jack == "end") t.jack_from_end = true;
            else if (jack != "start") throw py::value_error("tendon jack 'start', 'end' or 'both'");
            t.losses.mu = mu;
            t.losses.k = k;
            t.losses.slip = slip;
            t.losses.Ep = Ep;
            t.losses.Ap = Ap;
            t.losses.fpk = fpk;
            t.losses.rho1000 = rho1000;
            t.losses.relaxation_class = relaxation_class;
            t.s_begin = s_begin;
            t.s_end = s_end;
            for (const auto& [s0, s1] : free_runs)
              if (!(s1 > s0)) throw py::value_error("tendon free_runs: intervals (s0, s1) with s1 > s0");
            for (std::size_t i = 1; i < force_profile.size(); ++i)
              if (!(force_profile[i].first > force_profile[i - 1].first))
                throw py::value_error("tendon force_profile: points (s, P) with s strictly increasing");
            t.free_runs = std::move(free_runs);
            t.force_profile = std::move(force_profile);
            self.m.tendons.push_back(t);
            self.m.invalidate_index();
            return t.id;
          },
          py::arg("bars"), py::arg("P0"), py::arg("profile") = std::vector<std::array<double, 3>>{}, py::arg("jack") = "start", py::arg("mu") = 0.19,
          py::arg("k") = 0.01, py::arg("slip") = 0.006, py::arg("Ep") = 195e9, py::arg("Ap") = 0.0,
          py::arg("fpk") = 0.0, py::arg("rho1000") = 2.5, py::arg("relaxation_class") = 2,
          py::arg("s_begin") = 0.0, py::arg("s_end") = -1.0, py::arg("name") = "", py::arg("segments") = py::none(),
          py::arg("exact_angle") = false, py::arg("phi_abs") = 0.0, py::arg("non_simul") = 1,
          py::arg("relaxation") = "ec2", py::arg("relax_teb_days") = 0.0, py::arg("thermal_loss") = 0.0,
          py::arg("free_runs") = std::vector<std::pair<double, double>>{},
          py::arg("force_profile") = std::vector<std::pair<double, double>>{},
          "Add a prestressing tendon through a chain of bars: jack force P0 [N] at the active anchor(s) ('start', "
          "'end', 'both'), eccentricity profile [(s, e_y, e_z), ...] piecewise linear along the chain (s from the "
          "chain origin [m], eccentricities in each bar's section frame [m]), friction mu [1/rad], wobble k "
          "[rad/m], draw-in slip [m], tendon Ep, Ap [m^2] (0: force-only), fpk [Pa] and relaxation (rho1000 [%], "
          "class 1/2/3) for the deferred losses. Returns the tendon id; tensioned by the phases. segments = "
          "[dict(s0, s1, cy=(c0, c1, c2, c3), cz=(...))] replaces profile by contiguous cubic pieces e(s) = c0 + c1 ds "
          "+ c2 ds^2 + c3 ds^3, ds = s - s0 (exact curvature: parabolas, tangents); exact_angle: N = P cos(slope) "
          "instead of the small-slope N = P; phi_abs [1/m]: line friction exp(-phi_abs s) independent of mu; "
          "non_simul = n units tensioned one after the other ((n - 1) / 2n of their own elastic loss); relaxation "
          "'ec2' (EN 1992-1-1 §3.3.2) or 'bpel' (BPEL 91 annex 2, equivalent time); relax_teb_days, thermal_loss "
          "[Pa]: heat-cured pretension (EN 1992-1-1 §10.3.2). free_runs = [(s0, s1)]: chain intervals where the "
          "tendon runs free of its duct (external tendon), the wobble term k not acting there; force_profile = "
          "[(s, P)]: the force along the chain given point by point (piecewise linear, constant outside), which "
          "replaces friction, wobble and draw-in.")
      .def(
          "phased",
          [](const PyModel& self, const py::list& phases, const std::string& creep_engine, int stations,
             bool exact_losses, int substeps, bool creep_along_bar, bool origin_ledgers,
             std::optional<py::function> modulus, bool keep_structures) {
            PhasePlan plan;
            std::size_t i = 0;
            for (const py::handle h : phases) plan.phases.push_back(make_phase(self, h.cast<py::dict>(), i++));
            PhasingOptions o;
            if (creep_engine == "exact") o.creep_engine = CreepEngine::Exact;
            else if (creep_engine == "kelvin") o.creep_engine = CreepEngine::Kelvin;
            else if (creep_engine == "aaem") o.creep_engine = CreepEngine::Aaem;
            else if (creep_engine == "incremental" || creep_engine == "st1") o.creep_engine = CreepEngine::St1;
            else throw py::value_error("creep_engine 'exact', 'kelvin', 'aaem' or 'incremental'");
            o.stations_per_bar = stations;
            o.exact_losses = exact_losses;
            o.creep_substeps_per_decade = substeps;
            o.creep_along_bar = creep_along_bar;
            o.origin_ledgers = origin_ledgers;
            o.keep_structures = keep_structures;
            if (modulus) {
              py::function f = *modulus;
              o.modulus_at = [f](Id bar, double age) {
                py::gil_scoped_acquire g;
                const double E = f(bar, age).cast<double>();
                if (!(E > 0) || !std::isfinite(E)) throw Error("phased: modulus(bar, age) must return a finite E > 0");
                return E;
              };
            }
            ModelPtr model = self.solver_model();
            auto r = engine([&] {
              py::gil_scoped_release nogil;
              return std::make_shared<const PhasedResults>(solve_phased(*model, plan, o));
            });
            return PyPhased{model, r};
          },
          py::arg("phases"), py::arg("creep_engine") = "exact", py::arg("stations") = 11,
          py::arg("exact_losses") = true, py::arg("creep_substeps_per_decade") = 4,
          py::arg("creep_along_bar") = false, py::arg("origin_ledgers") = false, py::arg("modulus") = py::none(),
          py::arg("keep_structures") = false,
          "Staged construction: a list of phase dicts, in date order. Keys: name, time [days], activate / "
          "deactivate [bars], add_supports [dict(node, dx.., rz)], remove_supports [nodes], cases [load cases "
          "applied and kept], tension [dict(tendon, fraction=1, mode='jack'|'bed'|'release', unbonded=False, "
          "group=0)], inject [tendons], inject_partial [(tendon, s0, s1)], activate_age {bar: age}, regime "
          "('linear' | 'second_order'), imperfection [(node, (6 values))], modify_sections [(bar, section)], "
          "place_hinges / close_hinges [dict(bar, end=1, dofs='rz')], cable_actions [dict(cable, kind, value)]. "
          "Creep and shrinkage act between the phase dates. modulus = f(bar, age_days) -> E [Pa]: the instantaneous "
          "modulus of each increment at the age of the bar's concrete on the phase date (bars with a concrete law; "
          "default: the material's E). creep_engine 'exact' (default), 'kelvin', 'aaem' or 'incremental'. "
          "keep_structures=True records each phase's active substructure (PhaseResult.tangent()).");

  py::class_<PyPhased>(mod, "Phased", "Results of a staged-construction analysis, one entry per phase.")
      .def("__len__", [](const PyPhased& s) { return s.r->phases.size(); })
      .def(
          "phase",
          [](const PyPhased& s, std::size_t i) {
            if (i >= s.r->phases.size()) throw py::key_error("no phase " + std::to_string(i));
            return PyPhase{s.model, s.r, i};
          },
          py::arg("index"))
      .def("names", [](const PyPhased& s) {
        std::vector<std::string> v;
        for (const PhaseResult& p : s.r->phases) v.push_back(p.name);
        return v;
      });

  py::class_<PyPhase>(mod, "PhaseResult", "The state at the end of one phase (cumulated) and its increment.")
      .def_property_readonly("name", [](const PyPhase& s) { return s.get().name; })
      .def_property_readonly("time", [](const PyPhase& s) { return s.get().time; })
      .def("displacement", [](const PyPhase& s, Id node) { return dof6(s.get().cumulated.node_disp.at(node_pos(*s.model, node))); })
      .def("increment_displacement",
           [](const PyPhase& s, Id node) { return dof6(s.get().incremental.node_disp.at(node_pos(*s.model, node))); })
      .def("bar_forces", [](const PyPhase& s, Id bar) { return stations_of(bar_of(s.get().cumulated, bar)); },
           "Stations (x, N, Vy, Vz, T, My, Mz) of the cumulated state.")
      .def("increment_bar_forces", [](const PyPhase& s, Id bar) { return stations_of(bar_of(s.get().incremental, bar)); })
      .def("reactions", [](const PyPhase& s) { return reactions_of(s.get().cumulated); })
      .def("tendons",
           [](const PyPhase& s) {
             std::vector<std::tuple<Id, double, double, double>> out;
             for (const TendonPhaseForce& t : s.get().tendons) out.emplace_back(t.tendon, t.jacked_to, t.loss, t.effective);
             return out;
           },
           "Per tendon: (id, jacked force, cumulated loss, effective force) [N].")
      .def("cables",
           [](const PyPhase& s) {
             std::vector<std::tuple<Id, bool, double, double, double, double>> out;
             for (const PhaseResult::CableState& c : s.get().cables) out.emplace_back(c.bar, c.active, c.L0, c.H, c.T1, c.T2);
             return out;
           },
           "Catenary cables: (solver bar, active, L0, H, T1, T2).")
      .def("__repr__", [](const PyPhase& s) {
        return "<keta.PhaseResult '" + s.get().name + "' t = " + std::to_string(s.get().time) + " d>";
      });
}
