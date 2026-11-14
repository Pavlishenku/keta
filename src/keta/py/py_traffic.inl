// keta Python bindings — influence lines, traffic load models, fatigue,
// railway actions, transverse distribution, intensity laws, influence
// surfaces on node clouds. Included by keta_py.cpp inside its
// anonymous namespace.
//
//   il = m.influence_line(bars, dict(kind="bar_force", bar=b, x=2.5, comp=5))
//   keta.convoy_extreme(il, keta.flm3_vehicle())          # (max, min, r_max, r_min)
//   h = keta.convoy_history(il, keta.flm2_lorry(1))       # (r, effect)
//   keta.rainflow(h[1])                                   # [(range, mean, count)]
// A convoy is a dict(axles=[(offset, load)], segments=[(start, length, total)]).

ResponseKind response_of(const py::dict& response, const char* who) {
  const std::string kind = response.contains("kind") ? response["kind"].cast<std::string>() : "displacement";
  if (kind == "displacement") return ResponseKind::displacement(response["node"].cast<Id>(), response["comp"].cast<int>());
  if (kind == "reaction") return ResponseKind::reaction(response["node"].cast<Id>(), response["comp"].cast<int>());
  if (kind == "bar_force")
    return ResponseKind::bar_force(response["bar"].cast<Id>(), response["x"].cast<double>(), response["comp"].cast<int>());
  if (kind == "bar_displacement")
    return ResponseKind::bar_displacement(response["bar"].cast<Id>(), response["x"].cast<double>(),
                                          response["comp"].cast<int>());
  throw py::value_error(std::string(who) +
                        ": response kind 'displacement', 'reaction', 'bar_force' or 'bar_displacement'");
}

ConvoyPattern convoy_of(const py::dict& c) {
  ConvoyPattern p;
  if (c.contains("axles"))
    for (const py::handle h : c["axles"]) p.axles.push_back(h.cast<std::pair<double, double>>());
  if (c.contains("segments"))
    for (const py::handle h : c["segments"]) p.segments.push_back(h.cast<std::array<double, 3>>());
  for (const auto& [k, v] : c)
    if (std::string(py::str(k)) != "axles" && std::string(py::str(k)) != "segments")
      throw py::value_error("convoy: keys 'axles' [(offset, load)] and 'segments' [(start, length, total)]");
  return p;
}
py::dict convoy_dict(const ConvoyPattern& p) {
  py::dict d;
  d["axles"] = p.axles;
  d["segments"] = p.segments;
  return d;
}
py::tuple extremes_tuple(const Extremes& e) { return py::make_tuple(e.max, e.min, e.r_max, e.r_min); }

void bind_traffic(py::module_& mod) {
  auto cls = py::reinterpret_borrow<py::class_<PyModel>>(mod.attr("Model"));

  py::class_<InfluenceLine>(mod, "InfluenceLine",
                            "An influence line along a path of bars: abscissae s [m], ordinates eta (response per "
                            "unit load), cubic between samples when the slopes are known.")
      .def_readonly("s", &InfluenceLine::s)
      .def_readonly("eta", &InfluenceLine::eta)
      .def("at", &InfluenceLine::at, py::arg("s"), "Ordinate at s (0 off the path).")
      .def("integral", &InfluenceLine::integral, py::arg("a"), py::arg("b"))
      .def("integral_pos", &InfluenceLine::integral_pos, "Area of the positive zones.")
      .def("integral_neg", &InfluenceLine::integral_neg, "Area of the negative zones.")
      .def("point_extreme", &InfluenceLine::point_extreme, "(max, min) under a single unit load.")
      .def("zones",
           [](const InfluenceLine& il, int sign) {
             std::vector<std::tuple<double, double, double>> out;
             for (const LoadedZone& z : engine([&] { return sign_zones(il, sign); })) out.emplace_back(z.a, z.b, z.area);
             return out;
           },
           py::arg("sign") = 1, "The zones of one sign: [(s_start, s_end, area)].")
      .def_static(
          "from_samples",
          [](std::vector<double> s, std::vector<double> eta) {
            if (s.size() != eta.size() || s.size() < 2) throw py::value_error("influence line: matching s, eta (>= 2)");
            InfluenceLine il;
            il.s = std::move(s);
            il.eta = std::move(eta);
            return il;
          },
          py::arg("s"), py::arg("eta"), "A piecewise-linear line from samples (hand lines, tests).");

  cls.def(
      "influence_line",
      [](const PyModel& self, const std::vector<Id>& bars, const py::dict& response, std::array<double, 3> dir,
         int per_bar) {
        ModelPtr model = self.solver_model();
        const ResponseKind rk = response_of(response, "influence_line");
        return engine([&] {
          const LoadPath path = build_path(*model, bars, per_bar);
          return influence_line(*model, path, Vec3(dir[0], dir[1], dir[2]), rk);
        });
      },
      py::arg("bars"), py::arg("response"), py::arg("dir") = std::array<double, 3>{0, 0, -1}, py::arg("per_bar") = 11,
      "Influence line of one response (dict: kind 'displacement' / 'reaction' (node, comp), 'bar_force' / "
      "'bar_displacement' (bar, x, comp) — comp 0..5 = N, Vy, Vz, T, My, Mz for a bar force) for a unit load of "
      "direction `dir` travelling along the chain of `bars` (per_bar stations each). One factorization, the "
      "adjoint sweep.");
  cls.def(
      "influence_surface_cloud",
      [](const PyModel& self, const std::vector<Id>& nodes, const py::dict& response, std::array<double, 3> dir) {
        ModelPtr model = self.solver_model();
        const ResponseKind rk = response_of(response, "influence_surface_cloud");
        return engine([&] { return influence_surface_scattered(*model, nodes, Vec3(dir[0], dir[1], dir[2]), rk); });
      },
      py::arg("nodes"), py::arg("response"), py::arg("dir") = std::array<double, 3>{0, 0, -1},
      "Influence surface on ANY set of deck nodes (skew, fan-shaped): Delaunay triangulation of their plan "
      "positions, linear on each triangle, zero outside the hull.");

  py::class_<ScatteredSurface>(mod, "CloudSurface", "Influence surface on a node cloud (linear per triangle).")
      .def("at", &ScatteredSurface::at, py::arg("x"), py::arg("y"))
      .def("integral", &ScatteredSurface::integral)
      .def("integral_pos", &ScatteredSurface::integral_pos)
      .def("integral_neg", &ScatteredSurface::integral_neg)
      .def("wheels_value",
           [](const ScatteredSurface& s, const std::vector<std::array<double, 3>>& wheels, double rx, double ry) {
             WheelSet w;
             w.wheels = wheels;
             return wheel_set_value(s, w, rx, ry);
           },
           py::arg("wheels"), py::arg("rx"), py::arg("ry"), "Effect of wheels [(dx, dy, P)] at reference (rx, ry).")
      .def("wheels_extreme",
           [](const ScatteredSurface& s, const std::vector<std::array<double, 3>>& wheels, int sx, int sy) {
             WheelSet w;
             w.wheels = wheels;
             const SurfaceExtremes e = engine([&] { return wheel_set_extreme(s, w, sx, sy); });
             return py::make_tuple(e.max, e.min, py::make_tuple(e.rx_max, e.ry_max), py::make_tuple(e.rx_min, e.ry_min));
           },
           py::arg("wheels"), py::arg("steps_x") = 60, py::arg("steps_y") = 60,
           "(max, min, (x, y) of max, (x, y) of min) of the wheel set swept over the cloud.");

  // ---- moving loads on a line ----
  mod.def("convoy_extreme",
          [](const InfluenceLine& il, const py::dict& c, int steps) {
            return extremes_tuple(engine([&] { return convoy_extreme(il, convoy_of(c), steps); }));
          },
          py::arg("line"), py::arg("convoy"), py::arg("steps") = 500,
          "(max, min, r_max, r_min) of a convoy slid along an influence line (axles off the line count zero).");
  mod.def("udl_extreme", [](const InfluenceLine& il, double q) { return extremes_tuple(udl_extreme(il, q)); },
          py::arg("line"), py::arg("q"), "A UDL q free to load any zones: (q A+, q A-, 0, 0).");
  mod.def("convoy_history",
          [](const InfluenceLine& il, const py::dict& c, bool exact, int steps) {
            const EffectHistory h = engine([&] {
              return exact ? convoy_history_exact(il, convoy_of(c)) : convoy_history(il, convoy_of(c), steps);
            });
            return py::make_tuple(h.r, h.effect);
          },
          py::arg("line"), py::arg("convoy"), py::arg("exact") = true, py::arg("steps") = 500,
          "Effect history (r, E(r)) of a convoy crossing the line; exact=True samples every breakpoint (exact turning "
          "points for a piecewise-linear line), else a uniform sweep of steps + 1 positions.");
  mod.def("rainflow",
          [](const std::vector<double>& series) {
            std::vector<std::tuple<double, double, double>> out;
            for (const RainflowCycle& c : engine([&] { return rainflow(series); })) out.emplace_back(c.range, c.mean, c.count);
            return out;
          },
          py::arg("series"), "ASTM E1049 rainflow counting: [(range, mean, count)], count 1 or 0.5.");
  mod.def("turning_points", &turning_points, py::arg("series"));
  mod.def("multi_vehicle_extreme",
          [](const InfluenceLine& il, const py::list& vehicles, double udl, int steps) {
            std::vector<MobileVehicle> vs;
            for (const py::handle h : vehicles) {
              const py::dict v = h.cast<py::dict>();
              vs.push_back(MobileVehicle::from_margins(v["axles"].cast<std::vector<std::pair<double, double>>>(),
                                                       opt_get<double>(v, "occ_back", 0.0), opt_get<double>(v, "occ_front", 0.0),
                                                       opt_get<double>(v, "rem_back", 0.0), opt_get<double>(v, "rem_front", 0.0)));
            }
            const MultiVehicleResult r = engine([&] { return multi_vehicle_extreme(il, vs, udl, steps); });
            return py::make_tuple(r.extremes.max, r.extremes.min, r.max_r, r.min_r);
          },
          py::arg("line"), py::arg("vehicles"), py::arg("udl") = 0.0, py::arg("steps") = 200,
          "Several vehicles (dict axles, occ_back, occ_front, rem_back, rem_front) with disjoint occupancies and a "
          "concomitant UDL removed under them: the EXACT optimum on the position grid (max, min, positions at max, "
          "positions at min).");
  mod.def("intensity_udl_extreme",
          [](const InfluenceLine& il, py::function law, int sign) {
            const IntensityLaw L = IntensityLaw::custom([law](double l) {
              py::gil_scoped_acquire g;
              return law(l).cast<double>();
            });
            const IntensityExtreme e = engine([&] { return intensity_udl_extreme(il, L, sign); });
            return py::make_tuple(e.value, e.loaded_length, e.zones);
          },
          py::arg("line"), py::arg("density"), py::arg("sign") = 1,
          "UDL whose density q(L) depends on the loaded length: the exact worst subset of the zones of `sign` "
          "(value, loaded length, zone indices).");
  mod.def("courbon_coefficients", &courbon_coefficients, py::arg("girder_y"), py::arg("weights") = std::vector<double>{},
          py::arg("y"), "Courbon's transverse distribution: share of a unit load at y on each girder.");
  mod.def("separable_surface", &separable_surface, py::arg("line"), py::arg("ky"), py::arg("k"),
          "Influence surface IL(x) K(y) on the grid (line samples) x (ky).");

  // ---- road load models ----
  auto fac = [](const std::string& c) {
    if (c == "class1") return Lm1Factors::class1();
    if (c == "class2") return Lm1Factors::class2();
    throw py::value_error("factors 'class1' or 'class2'");
  };
  mod.def("en_lane_division",
          [](double w) {
            const LaneDivision d = en_lane_division(w);
            return py::make_tuple(d.n_lanes, d.lane_width, d.remaining);
          },
          py::arg("width"), "EN 1991-2 notional lanes: (n, lane width, remaining width).");
  mod.def("lm1_tandem", [fac](int rank, const std::string& f) { return convoy_dict(lm1_tandem(rank, fac(f))); },
          py::arg("rank"), py::arg("factors") = "class1");
  mod.def("lm1_udl", [fac](int rank, double w, const std::string& f) { return lm1_udl(rank, w, fac(f)); },
          py::arg("rank"), py::arg("lane_width"), py::arg("factors") = "class1", "UDL of a lane [N/m].");
  mod.def("lm1_lane_extreme",
          [fac](const InfluenceLine& il, int rank, double w, const std::string& f, int steps) {
            return extremes_tuple(engine([&] { return lm1_lane_extreme(il, rank, w, fac(f), steps); }));
          },
          py::arg("line"), py::arg("rank"), py::arg("lane_width"), py::arg("factors") = "class1", py::arg("steps") = 500);
  mod.def("lm2_axle", [](double beta, bool fr) { return convoy_dict(lm2_axle(beta, fr)); }, py::arg("beta_q") = 1.0,
          py::arg("french_na") = false);
  mod.def("lm4_crowd", &lm4_crowd, py::arg("lane_width"));
  mod.def("fasc61_a_of_l", &fasc61_a_of_l, py::arg("loaded_length"), "Fascicule 61 A(l) [N/m^2].");
  mod.def("fasc61_a1", &fasc61_a1, py::arg("bridge_class"), py::arg("loaded_lanes"));
  mod.def("fasc61_a2", &fasc61_a2, py::arg("bridge_class"), py::arg("lane_width"));
  mod.def("fasc61_a_lane_load", &fasc61_a_lane_load, py::arg("loaded_length"), py::arg("bridge_class"), py::arg("lane_width"));
  mod.def("fasc61_bc", [] { return convoy_dict(fasc61_bc()); });
  mod.def("fasc61_bt", [] { return convoy_dict(fasc61_bt()); });
  mod.def("fasc61_mc120", [] { return convoy_dict(fasc61_mc120()); });
  mod.def("fasc61_delta", &fasc61_delta, py::arg("element_length"), py::arg("permanent_load"), py::arg("max_axle_load"));
  // fatigue (EN 1991-2 §4.6)
  mod.def("flm1_tandem", [](int lane) { return convoy_dict(flm1_tandem(lane)); }, py::arg("lane"));
  mod.def("flm1_udl", &flm1_udl, py::arg("lane"), py::arg("lane_width"));
  mod.def("flm2_lorry", [](int i) { return convoy_dict(engine([&] { return flm2_lorry(i); })); }, py::arg("i"));
  mod.def("flm3_vehicle", [] { return convoy_dict(flm3_vehicle()); });
  mod.def("flm3_double", [](double gap) { return convoy_dict(engine([&] { return flm3_double(gap); })); }, py::arg("gap") = 40.0);
  mod.def("flm4_lorry", [](int i) { return convoy_dict(engine([&] { return flm4_lorry(i); })); }, py::arg("i"));
  mod.def("flm4_share",
          [](int i, const std::string& t) {
            const FatigueTraffic ft = t == "long" ? FatigueTraffic::LongDistance
                                      : t == "medium" ? FatigueTraffic::MediumDistance
                                      : t == "local" ? FatigueTraffic::Local
                                                     : throw py::value_error("traffic 'long', 'medium' or 'local'");
            return engine([&] { return flm4_share(i, ft); });
          },
          py::arg("i"), py::arg("traffic"), "Share [%] of FLM4 lorry i in the traffic type.");

  // ---- railway ----
  mod.def("rail_lm71", [](double a) { return convoy_dict(rail_lm71(a)); }, py::arg("alpha") = 1.0);
  mod.def("rail_sw0", [](double a) { return convoy_dict(rail_sw0(a)); }, py::arg("alpha") = 1.0);
  mod.def("rail_sw2", [] { return convoy_dict(rail_sw2()); });
  mod.def(
      "rail_hslm_b", [](int n_axles, double spacing) { return convoy_dict(engine([&] { return rail_hslm_b(n_axles, spacing); })); },
      py::arg("n_axles"), py::arg("spacing"),
      "HSLM-B of EN 1991-2 6.4.6.1.1(5): n_axles point forces of 170 kN at the uniform spacing [m], as a convoy "
      "dict (axles at 0, d, 2 d, ...). The number of forces and the spacing for a span are read by the caller on "
      "Figures 6.12 and 6.13 of the standard. For Model.rail_speed_sweep, give the axles with their offsets "
      "negated (the train behind its first axle).");
  mod.def("rail_unloaded_train", &rail_unloaded_train);
  mod.def("rail_phi2", &rail_phi2, py::arg("L_phi"));
  mod.def("rail_phi3", &rail_phi3, py::arg("L_phi"));
  auto model_of = [](const std::string& s) {
    if (s == "LM71") return RailModel::LM71;
    if (s == "SW0") return RailModel::SW0;
    if (s == "SW2") return RailModel::SW2;
    if (s == "HSLM") return RailModel::HSLM;
    throw py::value_error("rail model 'LM71', 'SW0', 'SW2' or 'HSLM'");
  };
  mod.def("rail_traction", &rail_traction, py::arg("loaded_length"), py::arg("alpha") = 1.0);
  mod.def("rail_braking",
          [model_of](double L, const std::string& m, double a) { return rail_braking(L, model_of(m), a); },
          py::arg("loaded_length"), py::arg("model") = "LM71", py::arg("alpha") = 1.0);
  mod.def("rail_nosing", &rail_nosing, py::arg("alpha") = 1.0);
  mod.def("rail_centrifugal_factor", &rail_centrifugal_factor, py::arg("v_kmh"), py::arg("L_f"));
  mod.def("rail_centrifugal_force", &rail_centrifugal_force, py::arg("v_kmh"), py::arg("radius"), py::arg("q_vertical"),
          py::arg("L_f"));
  mod.def("eccentric_rail",
          [](double total, double ratio, double gauge) {
            const EccentricRail e = engine([&] { return eccentric_rail(total, ratio, gauge); });
            return py::make_tuple(e.q1, e.q2, e.e());
          },
          py::arg("total"), py::arg("ratio") = 1.25, py::arg("gauge") = 1.5, "(q1, q2, eccentricity).");
  mod.def("ballast_patches",
          [](double Q, double depth, const std::string& sleeper, bool three, double spacing) {
            BallastOptions o;
            o.depth = depth;
            o.sleeper = sleeper == "monobloc" ? Sleeper::Monobloc : Sleeper::TwinBlock;
            o.three_sleepers = three;
            o.spacing = spacing;
            std::vector<std::tuple<double, double, double, double, double>> out;
            for (const LoadPatch& p : engine([&] { return ballast_patches(Q, o); })) out.emplace_back(p.dx, p.dy, p.load, p.ax, p.ay);
            return out;
          },
          py::arg("axle_load"), py::arg("depth") = 0.3, py::arg("sleeper") = "monobloc", py::arg("three_sleepers") = true,
          py::arg("spacing") = 0.6, "Footprints [(dx, dy, load, ax, ay)] of an axle through the ballast.");
  mod.def("annex_d_train",
          [](int n) {
            const ServiceTrain t = engine([&] { return annex_d_train(n); });
            py::dict d = convoy_dict(t.pattern);
            d["number"] = t.number;
            d["mass_t"] = t.mass_t;
            d["v_max_kmh"] = t.v_max_kmh;
            return d;
          },
          py::arg("number"), "EN 1991-2 Annex D service train 1..12: dict(axles, segments, number, mass_t, v_max_kmh).");
  mod.def("annex_d_mix",
          [](const std::string& s) {
            const RailFatigueMix m = s == "standard" ? RailFatigueMix::Standard
                                    : s == "heavy" ? RailFatigueMix::Heavy250
                                    : s == "light" ? RailFatigueMix::Light
                                                   : throw py::value_error("mix 'standard', 'heavy' or 'light'");
            std::vector<std::tuple<int, int, double>> out;
            for (const RailMixEntry& e : annex_d_mix(m)) out.emplace_back(e.train, e.per_day, e.tonnage_mt);
            return out;
          },
          py::arg("mix"), "[(train, trains per day, 1e6 t / year)].");
  mod.def("rail_n0_lower", &rail_n0_lower, py::arg("L_phi"));
  mod.def("rail_fatigue_dynamic_factor", &rail_fatigue_dynamic_factor, py::arg("v_ms"), py::arg("L_phi"), py::arg("n0"));
}
