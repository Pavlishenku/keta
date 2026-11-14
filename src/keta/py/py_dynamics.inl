// keta Python bindings — linear time histories (direct Newmark / HHT or
// modal superposition) under load-case time functions, a ground
// acceleration, a harmonic force or a moving convoy; floor response spectra;
// strain-energy modal damping; directional combinations. Included
// by keta_py.cpp inside its anonymous namespace.
//
//   h = m.time_history(dt=0.01, steps=500, cases=[(1, [0, 1], [1, 1])],
//                      damping=0.05, observe=[(n, 2)])
//   h["u"][(n, 2)], h["t"]

DampingSpec damping_of(const py::object& o) {
  DampingSpec d;
  if (o.is_none()) {  // undamped
    d.kind = DampingSpec::Kind::Rayleigh;
    return d;
  }
  if (py::isinstance<py::float_>(o) || py::isinstance<py::int_>(o)) {
    d.kind = DampingSpec::Kind::Modal;
    d.zeta = o.cast<double>();
    return d;
  }
  const py::dict m = o.cast<py::dict>();
  if (m.contains("per_mode")) {
    d.kind = DampingSpec::Kind::PerMode;
    d.per_mode = m["per_mode"].cast<std::vector<double>>();
  } else if (m.contains("rayleigh")) {
    d.kind = DampingSpec::Kind::Rayleigh;
    std::tie(d.alpha, d.beta) = m["rayleigh"].cast<std::pair<double, double>>();
  } else {
    throw py::value_error("damping: a float (modal ratio), dict(per_mode=[...]) or dict(rayleigh=(alpha, beta))");
  }
  return d;
}

IntegrationScheme scheme_of(const std::string& s, double hht_alpha) {
  if (s == "average") return IntegrationScheme::newmark_average();
  if (s == "linear") return IntegrationScheme::newmark_linear();
  if (s == "hht") return IntegrationScheme::hht(hht_alpha);
  throw py::value_error("scheme 'average', 'linear' or 'hht'");
}

std::function<double(double)> table_fn(std::vector<double> ts, std::vector<double> vs) {
  if (ts.size() != vs.size() || ts.empty()) throw py::value_error("time function: times and values of one length");
  return [ts, vs](double t) {
    if (t <= ts.front()) return vs.front();
    if (t >= ts.back()) return vs.back();
    const std::size_t k = static_cast<std::size_t>(std::upper_bound(ts.begin(), ts.end(), t) - ts.begin());
    const double a = (t - ts[k - 1]) / (ts[k] - ts[k - 1]);
    return vs[k - 1] + a * (vs[k] - vs[k - 1]);
  };
}

PedestrianDirection pedestrian_direction_of(const std::string& s) {
  if (s == "vertical") return PedestrianDirection::Vertical;
  if (s == "longitudinal") return PedestrianDirection::Longitudinal;
  if (s == "lateral") return PedestrianDirection::Lateral;
  throw py::value_error("direction 'vertical', 'longitudinal' or 'lateral'");
}

FootbridgeClass footbridge_class_of(const py::object& o) {
  if (py::isinstance<py::int_>(o)) {
    const int c = o.cast<int>();
    if (c >= 1 && c <= 4) return static_cast<FootbridgeClass>(c);
  } else if (py::isinstance<py::str>(o)) {
    const std::string c = o.cast<std::string>();
    if (c == "I") return FootbridgeClass::I;
    if (c == "II") return FootbridgeClass::II;
    if (c == "III") return FootbridgeClass::III;
    if (c == "IV") return FootbridgeClass::IV;
  }
  throw py::value_error("footbridge class 1..4 or 'I', 'II', 'III', 'IV'");
}

struct PyTangent {
  ModelPtr source;  // the model the state came from (consumers check its fingerprint)
  std::shared_ptr<const TangentState> t;
};

PyTangent make_tangent(ModelPtr source, TangentState t) {
  return PyTangent{std::move(source), std::make_shared<const TangentState>(std::move(t))};
}

void bind_dynamics(py::module_& mod) {
  auto M = py::reinterpret_borrow<py::class_<PyModel>>(mod.attr("Model"));

  M.def(
      "time_history",
      [](PyModel& self, double dt, int steps, const std::vector<std::tuple<Id, std::vector<double>, std::vector<double>>>& cases,
         std::optional<std::pair<int, std::vector<double>>> ground,
         std::optional<std::tuple<Id, int, double, double>> harmonic, const py::object& damping,
         const std::string& scheme, double hht_alpha, const py::object& modal, bool lumped,
         const std::vector<std::pair<Id, int>>& observe) {
        ModelPtr model = self.solver_model();
        const Model& md = *model;
        Assembled a = engine([&] { return assemble_elements(md, true); });
        // the load-case vectors, scaled by their time functions
        std::vector<std::pair<VecX, std::function<double(double)>>> parts;
        for (const auto& [cid, ts, vs] : cases) {
          const LoadCase* lc = nullptr;
          for (const LoadCase& c : md.cases)
            if (c.id == cid) lc = &c;
          if (!lc) throw py::key_error("no load case " + std::to_string(cid));
          CaseLoads cl = engine([&] { return lower_case(md, a, *lc); });
          parts.emplace_back(cl.rhs, table_fn(ts, vs));
        }
        std::function<VecX(int)> seis, harm;
        if (ground) {
          if (ground->first < 0 || ground->first > 2) throw py::value_error("ground direction 0, 1 or 2");
          seis = engine([&] { return seismic_force(md, ground->first, ground->second, lumped); });
        }
        if (harmonic) {
          const auto& [node, comp, amp, omega] = *harmonic;
          harm = engine([&] { return harmonic_force(md, node, comp, amp, omega, dt); });
        }
        const int n = a.dofs.n_free;
        auto force = [&, n](int step) {
          VecX f = VecX::Zero(n);
          for (const auto& [v, fn] : parts) f += v * fn(step * dt);
          if (seis) f += seis(step);
          if (harm) f += harm(step);
          return f;
        };
        const DampingSpec dmp = damping_of(damping);
        const IntegrationScheme sch = scheme_of(scheme, hht_alpha);
        TimeHistoryResults r = engine([&] {
          if (!modal.is_none())
            return time_history_modal(md, *modal.cast<const PyModal&>().r, steps, dt, force, dmp, sch);
          return time_integrate(md, steps, dt, force, dmp, sch, lumped);
        });
        py::dict out;
        std::vector<double> t;
        for (std::size_t i = 0; i < r.u.size(); ++i) t.push_back(static_cast<double>(i) * dt);
        out["t"] = t;
        out["energy"] = r.energy;
        py::dict u, v, acc;
        for (const auto& [node, comp] : observe) {
          if (comp < 0 || comp > 5) throw py::value_error("observed component 0..5");
          const int ix = a.dofs.free_of(node_pos(md, node), comp);
          std::vector<double> hu(r.u.size(), 0.0), hv(r.u.size(), 0.0), ha(r.u.size(), 0.0);
          if (ix >= 0)
            for (std::size_t i = 0; i < r.u.size(); ++i) {
              hu[i] = r.u[i](ix);
              if (i < r.v.size()) hv[i] = r.v[i](ix);
              if (i < r.a.size()) ha[i] = r.a[i](ix);
            }
          const py::tuple key = py::make_tuple(node, comp);
          u[key] = hu;
          v[key] = hv;
          acc[key] = ha;
        }
        out["u"] = u;
        out["v"] = v;
        out["a"] = acc;
        return out;
      },
      py::kw_only(), py::arg("dt"), py::arg("steps"),
      py::arg("cases") = std::vector<std::tuple<Id, std::vector<double>, std::vector<double>>>{},
      py::arg("ground") = py::none(), py::arg("harmonic") = py::none(), py::arg("damping") = py::none(),
      py::arg("scheme") = "average", py::arg("hht_alpha") = 0.0, py::arg("modal") = py::none(),
      py::arg("lumped") = false, py::arg("observe") = std::vector<std::pair<Id, int>>{},
      "Linear time history, direct (Newmark / HHT on M, C, K) or by modal superposition (modal = a "
      "Model.modal(...) result). Excitations, summed: cases = [(case id, times, values)] (each case's load "
      "vector times its piecewise-linear time function), ground = (direction 0|1|2, accelerations per step) "
      "(relative displacements), harmonic = (node, comp, amplitude, omega). damping: None (undamped), a modal "
      "ratio or dict(per_mode=[...]) (modal superposition only), dict(rayleigh=(alpha, beta)) (both). Returns dict(t, energy, u / v / a = {(node, comp): "
      "history}) for the observed DOFs (0 on a blocked DOF).");

  M.def(
      "floor_spectrum",
      [](PyModel& self, const PyModal& modal, Id node, int comp, int dir, const std::vector<double>& accel, double dt,
         const std::vector<double>& periods, double zeta_floor, const py::object& damping) {
        const FloorSpectrum f = engine([&] {
          return floor_response_spectrum(*modal.model, *modal.r, node, comp, dir, accel, dt, damping_of(damping),
                                         periods, zeta_floor);
        });
        (void)self;
        py::dict d;
        d["accel"] = f.accel;
        d["periods"] = f.periods;
        d["sa"] = f.sa;
        return d;
      },
      py::arg("modal"), py::arg("node"), py::arg("comp"), py::arg("dir"), py::arg("accel"), py::arg("dt"),
      py::arg("periods"), py::arg("zeta_floor") = 0.05, py::arg("damping") = py::none(),
      "Floor response spectrum: the absolute acceleration history of (node, comp) under a ground acceleration "
      "along `dir` by modal superposition (exact steps), then its pseudo-acceleration spectrum at zeta_floor. "
      "Returns dict(accel, periods, sa).");

  M.def(
      "strain_energy_damping",
      [](PyModel& self, const PyModal& modal, const py::object& zeta, double zeta_links) {
        (void)self;
        std::map<Id, double> zeta_of_material;
        if (zeta.is_none())
          for (const Material& mat : modal.model->materials) zeta_of_material[mat.id] = mat.damping;
        else
          zeta_of_material = zeta.cast<std::map<Id, double>>();
        return engine([&] {
          return strain_energy_damping(*modal.model, *modal.r, zeta_of_material, zeta_links).per_mode;
        });
      },
      py::arg("modal"), py::arg("zeta") = py::none(), py::arg("zeta_links") = 0.0,
      "Composite modal damping zeta_i = sum_e zeta_e E_e,i / sum_e E_e,i (modal strain energies of bars, shells "
      "and elastic links); zeta = {material id: ratio}, or None for the damping ratio each material was given "
      "(add_material(damping=)). Returns the per-mode ratios (pass dict(per_mode=...)).");

  mod.def(
      "combine_directional",
      [](const std::vector<std::vector<Dof6>>& dirs, const std::string& rule) {
        std::vector<std::vector<Vec6>> v;
        for (const auto& field : dirs) {
          std::vector<Vec6> f;
          for (const Dof6& x : field) {
            Vec6 e;
            e << std::get<0>(x), std::get<1>(x), std::get<2>(x), std::get<3>(x), std::get<4>(x), std::get<5>(x);
            f.push_back(e);
          }
          v.push_back(f);
        }
        std::vector<Vec6> r;
        if (rule == "srss") r = engine([&] { return combine_directional_srss(v); });
        else if (rule == "30%") r = engine([&] { return combine_directional_30pct(v); });
        else throw py::value_error("directional rule 'srss' or '30%'");
        std::vector<Dof6> out;
        for (const Vec6& x : r) out.push_back(dof6(x));
        return out;
      },
      py::arg("fields"), py::arg("rule") = "srss",
      "Combines per-direction peak fields (lists of 6-tuples, one list per direction): 'srss' or the EN 1998-1 "
      "30 % rule (max over the leading direction of E_lead + 0.3 sum of the others).");

  // ---- accidental rupture during a time history ----
  M.def(
      "rupture_history",
      [](PyModel& self, double dt, int steps, int rupture_step, const std::vector<Id>& bars,
         const std::vector<std::tuple<Id, std::vector<double>, std::vector<double>>>& cases, const py::object& damping,
         const std::string& scheme, double hht_alpha, bool lumped, bool start_static,
         const std::vector<std::pair<Id, int>>& observe) {
        ModelPtr model = self.solver_model();
        const Model& md = *model;
        RuptureHistoryOptions o;
        o.dt = dt;
        o.steps = steps;
        o.rupture_step = rupture_step;
        o.bars = bars;
        o.damping = damping_of(damping);
        if (o.damping.kind != DampingSpec::Kind::Rayleigh)
          throw py::value_error("rupture_history: damping None or dict(rayleigh=(alpha, beta)) (direct integration)");
        o.scheme = scheme_of(scheme, hht_alpha);
        o.lumped_mass = lumped;
        o.start_static = start_static;
        std::vector<TimedCase> tc;
        for (const auto& [cid, ts, vs] : cases) tc.push_back({cid, table_fn(ts, vs)});
        const RuptureHistoryResults r = engine([&] { return rupture_time_history(md, tc, o); });
        py::dict out;
        std::vector<double> t;
        for (std::size_t i = 0; i < r.u.size(); ++i) t.push_back(static_cast<double>(i) * dt);
        out["t"] = t;
        out["rupture_step"] = r.rupture_step;
        py::dict u, v, acc, s0, s1;
        for (const auto& [node, comp] : observe) {
          if (comp < 0 || comp > 5) throw py::value_error("observed component 0..5");
          const std::size_t ni = node_pos(md, node);
          std::vector<double> hu, hv, ha;
          for (std::size_t i = 0; i < r.u.size(); ++i) {
            hu.push_back(r.u[i][ni](comp));
            hv.push_back(r.v[i][ni](comp));
            ha.push_back(r.a[i][ni](comp));
          }
          const py::tuple key = py::make_tuple(node, comp);
          u[key] = hu;
          v[key] = hv;
          acc[key] = ha;
          s0[key] = r.static_intact[ni](comp);
          s1[key] = r.static_severed[ni](comp);
        }
        out["u"] = u;
        out["v"] = v;
        out["a"] = acc;
        out["static_intact"] = s0;
        out["static_severed"] = s1;
        return out;
      },
      py::kw_only(), py::arg("dt"), py::arg("steps"), py::arg("rupture_step"), py::arg("bars"),
      py::arg("cases") = std::vector<std::tuple<Id, std::vector<double>, std::vector<double>>>{},
      py::arg("damping") = py::none(), py::arg("scheme") = "average", py::arg("hht_alpha") = 0.0,
      py::arg("lumped") = false, py::arg("start_static") = true,
      py::arg("observe") = std::vector<std::pair<Id, int>>{},
      "Linear time history in which `bars` vanish at t = rupture_step x dt (a cut stay, suspender or strut): "
      "the intact structure to the rupture, then the structure without them resumed from the state reached "
      "(displacements and velocities continuous, the acceleration jumps). cases = [(case id, times, values)]; "
      "start_static: t = 0 at the intact static equilibrium under the loads at t = 0. damping: None or "
      "dict(rayleigh=(alpha, beta)). Returns dict(t, rupture_step, u / v / a = {(node, comp): history}, "
      "static_intact / static_severed = {(node, comp): value}).");

  // ---- saved tangent states ----
  py::class_<PyTangent>(mod, "Tangent",
                        "The structure linearized about an equilibrium state (saved tangent): modal and buckling "
                        "analyses about that state.")
      .def_property_readonly("source",
                             [](const PyTangent& s) {
                               switch (s.t->source) {
                                 case TangentSource::StaticCase: return std::string("static");
                                 case TangentSource::CableState: return std::string("cable");
                                 case TangentSource::Corotational: return std::string("corotational");
                                 case TangentSource::PhasingStage: return std::string("phase");
                               }
                               return std::string();
                             })
      .def_property_readonly("label", [](const PyTangent& s) { return s.t->label; })
      .def_property_readonly("case_id", [](const PyTangent& s) { return s.t->case_id; })
      .def_property_readonly("converged", [](const PyTangent& s) { return s.t->converged; })
      .def_property_readonly("asymmetry", [](const PyTangent& s) { return s.t->asymmetry; })
      .def_property_readonly("bar_ids",
                             [](const PyTangent& s) {
                               std::vector<Id> ids;
                               for (const Bar& b : s.t->structure->bars) ids.push_back(b.id);
                               return ids;
                             })
      .def(
          "axial",
          [](const PyTangent& s, Id bar) {
            const int i = s.t->structure->bar_index(bar);
            if (i < 0) throw py::key_error("bar " + std::to_string(bar) + " is not in the linearized structure");
            return s.t->axial[static_cast<std::size_t>(i)];
          },
          py::arg("bar"), "Axial force N of a bar in the state [N] (tension +; 0 on catenary cables).")
      .def(
          "displacement",
          [](const PyTangent& s, Id node) {
            const int i = s.t->structure->node_index(node);
            return dof6(s.t->u[static_cast<std::size_t>(i)]);
          },
          py::arg("node"), "Displacement of a node in the state (6-tuple).")
      .def(
          "modal",
          [](const PyTangent& s, int n_modes, bool lumped) {
            ModalOptions o;
            o.n_modes = n_modes;
            o.lumped_mass = lumped;
            auto r = engine([&] { return std::make_shared<const ModalResults>(solve_modal_on_state(*s.source, *s.t, o)); });
            return PyModal{s.t->structure, r};
          },
          py::arg("n_modes") = 6, py::arg("lumped") = false,
          "Modal analysis on the saved tangent: the modes of the loaded (or staged) structure. Tension stiffens, "
          "compression softens; the shapes are per node of the linearized structure.")
      .def(
          "buckling",
          [](const PyTangent& s, const PyModel& model, Id cs, int n_modes) {
            ModelPtr m = model.solver_model();
            const LoadCase* lc = nullptr;
            for (const LoadCase& c : m->cases)
              if (c.id == cs) lc = &c;
            if (!lc) throw py::key_error("no load case " + std::to_string(cs));
            auto r = engine([&] {
              return std::make_shared<const BucklingResults>(solve_buckling_on_state(*m, *s.t, lc->loads, n_modes));
            });
            return PyBuckling{s.t->structure, r};
          },
          py::arg("model"), py::arg("case"), py::arg("n_modes") = 3,
          "Buckling of the ADDITIONAL loads of case `case` of `model` about the saved state: K_T + lambda "
          "Kg(N_case) singular (factors scale the additional case only; the state stays as it is). `model` is the "
          "model the state came from (its structure is checked by fingerprint; its cases may have grown).");
  M.def(
      "tangent",
      [](const PyModel& self, std::optional<Id> cs) {
        ModelPtr model = self.solver_model();
        const Id id = self.case_or_first(cs);
        return make_tangent(model, engine([&] { return tangent_of_static_case(*model, id); }));
      },
      py::arg("case") = py::none(),
      "Saved tangent about the linear static solution of a case: K + Kg(N) (bars and shell membranes).");
  py::reinterpret_borrow<py::class_<PyNonlinear>>(mod.attr("Nonlinear"))
      .def(
          "tangent",
          [](const PyNonlinear& s) {
            return make_tangent(s.model, engine([&] {
                                  return tangent_of_corotational(*s.model, s.case_id, *s.r, CorotRotationMeasure::Sin);
                                }));
          },
          "Saved tangent about the converged large-displacement state: the consistent corotational tangent "
          "(catenary cables included).");
  py::reinterpret_borrow<py::class_<PyCableStatic>>(mod.attr("CableStatic"))
      .def(
          "tangent",
          [](const PyCableStatic& s) {
            return make_tangent(s.model, engine([&] { return tangent_of_cable_state(*s.model, s.case_id, *s.r); }));
          },
          "Saved tangent about the cable-static state: frame K + Kg(N) plus the catenary pair tangents at the "
          "deformed anchors.");
  py::reinterpret_borrow<py::class_<PyPhase>>(mod.attr("PhaseResult"))
      .def(
          "tangent",
          [](const PyPhase& s) { return make_tangent(s.model, engine([&] { return tangent_of_phase(*s.model, s.get()); })); },
          "Saved tangent about the end of the phase: its active substructure (phased(keep_structures=True)) at "
          "the moduli of its date, K + Kg of the cumulated axial field.");

  // ---- footbridges under pedestrians (Setra / AFGC 2006) ----
  M.def(
      "pedestrian_load",
      [](PyModel& self, const PyModal& modal, int mode, double q, int comp, std::optional<std::vector<Id>> bars,
         Id cs, bool sign_cancelled) {
        PedestrianModeLoad pl;
        pl.q = q;
        pl.comp = comp;
        pl.sign_cancelled = sign_cancelled;
        if (bars) pl.bars = *bars;
        const LoadCase lc = engine([&] { return pedestrian_mode_load_case(*modal.model, *modal.r, mode, pl); });
        for (const Load& l : lc.loads) (void)self.bar(l.bar);  // the bars of the solved model exist here
        LoadCase& c = self.case_mut(cs);
        c.loads.insert(c.loads.end(), lc.loads.begin(), lc.loads.end());
        return lc.loads.size();
      },
      py::arg("modal"), py::arg("mode"), py::arg("q"), py::arg("comp") = 2, py::arg("bars") = py::none(),
      py::arg("case") = 1, py::arg("sign_cancelled") = false,
      "Adds to load case `case` the pedestrian load that follows mode `mode` (0-based) of `modal`: a uniform global "
      "load q sign(phi) [N/m] along global direction comp (0, 1, 2) on every bar (or on `bars`), the sign of the "
      "mode-shape component along the bar, each bar cut at the zeros of its shape (Setra 2006: every contribution "
      "unfavourable). sign_cancelled=True keeps q only "
      "where phi > 0. Drive it at the natural frequency with harmonic() or time_history(cases=...). Returns the number "
      "of bar loads added.");

  mod.def(
      "setra_equivalent_pedestrians",
      [](double n, double xi, bool very_dense) { return engine([&] { return setra_equivalent_pedestrians(n, xi, very_dense); }); },
      py::arg("n"), py::arg("xi"), py::arg("very_dense") = false,
      "Setra 2006 equivalent number of pedestrians in phase: 10.8 sqrt(xi n) (random frequencies and phases), "
      "1.85 sqrt(n) (very dense crowd, random phases).");
  mod.def(
      "setra_frequency_range",
      [](double f, const std::string& direction) {
        return engine([&] { return setra_frequency_range(f, pedestrian_direction_of(direction)); });
      },
      py::arg("f"), py::arg("direction") = "vertical",
      "Setra 2006 frequency range 1 (maximum risk of resonance) to 4 (negligible) of a natural frequency [Hz] "
      "(tables 2.3 / 2.4).");
  mod.def(
      "setra_comfort_range",
      [](double accel, const std::string& direction) {
        return engine([&] { return setra_comfort_range(accel, pedestrian_direction_of(direction)); });
      },
      py::arg("accel"), py::arg("direction") = "vertical",
      "Setra 2006 acceleration range 1 (maximum comfort) to 4 (unacceptable) of a peak acceleration [m/s2] "
      "(tables 2.1 / 2.2: 0.5 / 1 / 2.5 vertical, 0.15 / 0.3 / 0.8 horizontal).");
  mod.attr("SETRA_LOCK_IN_ACCEL") = kSetraLockInAccel;
  mod.def(
      "setra_psi",
      [](double f, const std::string& direction, int harmonic) {
        return engine([&] { return setra_psi(f, pedestrian_direction_of(direction), harmonic); });
      },
      py::arg("f"), py::arg("direction") = "vertical", py::arg("harmonic") = 1,
      "Setra 2006 reduction factor psi(f) (figure 2.3 for harmonic 1, figure 2.4 for harmonic 2).");
  mod.def(
      "setra_crowd_density",
      [](const py::object& cls, int load_case) {
        return engine([&] { return setra_crowd_density(footbridge_class_of(cls), load_case); });
      },
      py::arg("footbridge_class"), py::arg("load_case"),
      "Setra 2006 crowd density [pedestrians/m2] of load case 1, 2 or 3 for a class (1..4 or 'I'..'IV').");
  mod.def(
      "setra_load_cases",
      [](const py::object& cls, int range) {
        return engine([&] { return setra_load_cases(footbridge_class_of(cls), range); });
      },
      py::arg("footbridge_class"), py::arg("frequency_range"),
      "Setra 2006 load cases to run for a mode in a frequency range (table 2.5).");
  mod.def(
      "setra_surface_load",
      [](const py::object& cls, int load_case, const std::string& direction, double f, double xi, double area) {
        const SetraSurfaceLoad s = engine([&] {
          return setra_surface_load(footbridge_class_of(cls), load_case, pedestrian_direction_of(direction), f, xi, area);
        });
        py::dict d;
        d["density"] = s.density;
        d["pedestrians"] = s.pedestrians;
        d["equivalent"] = s.equivalent;
        d["force"] = s.force;
        d["psi"] = s.psi;
        d["amplitude"] = s.amplitude;
        d["frequency"] = s.load_frequency;
        return d;
      },
      py::arg("footbridge_class"), py::arg("load_case"), py::arg("direction"), py::arg("f"), py::arg("xi"),
      py::arg("area"),
      "Setra 2006 harmonic surface load of a load case: amplitude = d F0 n_eq / N psi(f) [N/m2] with N = d area. "
      "Returns dict(density, pedestrians, equivalent, force, psi, amplitude, frequency).");
  mod.def(
      "fascicule61_footbridge_load",
      [](double loaded_length) { return engine([&] { return st1_al_pieton(loaded_length) * 1000.0 * kGravity; }); },
      py::arg("loaded_length"),
      "Uniform crowd load of a footbridge reserved to pedestrians and cycles [N/m2], Fascicule 61 titre II art. "
      "13.2: a(l) = 200 + 15000 / (l + 50) kg/m2 (kgf, 9.80665 N), l the loaded length [m].");
  mod.def("setra_pedestrian_mass", &setra_pedestrian_mass, py::arg("density"),
          "Pedestrian mass to add to the deck [kg/m2]: 70 kg per pedestrian.");
  mod.def(
      "fourier_load_factors",
      [](double a0, const std::vector<double>& a_cos, const std::vector<double>& b_sin, double period, double dt,
         int nsteps) { return engine([&] { return fourier_load_factors(a0, a_cos, b_sin, period, dt, nsteps); }); },
      py::arg("a0"), py::arg("a_cos"), py::arg("b_sin"), py::arg("period"), py::arg("dt"), py::arg("nsteps"),
      "Periodic load factor a0 + sum_i [a_i cos(2 pi i t / T) + b_i sin(2 pi i t / T)] sampled at t = n dt, "
      "n = 0..nsteps (a_cos[0], b_sin[0]: the first harmonic).");
  mod.def(
      "high_frequency_filtered_damping",
      [](const std::vector<double>& freqs, const std::vector<double>& xi, double f_c, double d_c) {
        return engine([&] { return high_frequency_filtered_damping(freqs, xi, f_c, d_c); });
      },
      py::arg("freqs"), py::arg("xi"), py::arg("f_c"), py::arg("d_c"),
      "Per-mode damping ratios xi_i, plus d_c above the cut frequency f_c [Hz] (pass dict(per_mode=...)).");
}
