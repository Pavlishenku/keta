// keta Python bindings — the isolated elastic catenary, tension control
// (target -> unstrained length) and Ernst's equivalent linear stiffness.
// Included by keta_py.cpp inside its anonymous namespace.
//
//   c = keta.catenary(dx=150, dz=80, l0=170, ea=2e9, w=800)   # dict(H, V1, V2, T1, T2)
//   l0 = keta.catenary_l0(150, 80, ea=2e9, w=800, end_tension=2.5e6, end=1)
//   k = keta.ernst_stiffness(ea=2e9, l0=200, dx=200, dz=0, w=800, h=2e6)["k"]

void bind_cables(py::module_& mod) {
  mod.def(
      "catenary",
      [](double dx, double dz, double l0, double ea, double w) {
        const Catenary2D c = engine([&] { return solve_catenary(dx, dz, l0, ea, w); });
        py::dict d;
        d["H"] = c.H;
        d["V1"] = c.V1;
        d["V2"] = c.V2();
        d["T1"] = c.T1();
        d["T2"] = c.T2();
        return d;
      },
      py::arg("dx"), py::arg("dz"), py::arg("l0"), py::arg("ea"), py::arg("w"),
      "The elastic catenary (Irvine) of unstrained length l0, rigidity ea [N] and weight w [N per unstrained m] "
      "spanning the chord dx > 0 (horizontal), dz (vertical, end 2 minus end 1): dict(H, V1, V2, T1, T2), the "
      "horizontal tension, the vertical components at both ends (V2 = V1 + w l0) and the end tensions [N].");
  mod.def(
      "catenary_l0",
      [](double dx, double dz, double ea, double w, std::optional<double> h, std::optional<double> end_tension,
         int end) {
        if (h.has_value() == end_tension.has_value())
          throw py::value_error("catenary_l0: give exactly one of h= (horizontal tension) or end_tension=");
        return engine([&] {
          return h ? catenary_l0_for_tension(dx, dz, ea, w, *h) : catenary_l0_for_end_tension(dx, dz, ea, w, *end_tension, end);
        });
      },
      py::arg("dx"), py::arg("dz"), py::kw_only(), py::arg("ea"), py::arg("w"), py::arg("h") = py::none(),
      py::arg("end_tension") = py::none(), py::arg("end") = 2,
      "Tension control: the unstrained length of the elastic catenary on the chord (dx, dz) whose horizontal "
      "tension is h, or whose tension at end `end` (1 or 2) is end_tension — on the taut branch (a target below "
      "the tension minimum of the hanging cable or above a 4.5 % strain is refused).");
  mod.def("ernst_modulus", [](double E, double gamma, double l, double sigma) {
            return engine([&] { return ernst_modulus(E, gamma, l, sigma); });
          },
          py::arg("e"), py::arg("gamma"), py::arg("l"), py::arg("sigma"),
          "Ernst's equivalent modulus E / (1 + (gamma l)^2 E / (12 sigma^3)): specific weight gamma = w / A [N/m3], "
          "horizontal projection l [m], tension stress sigma [Pa].");
  mod.def(
      "ernst_stiffness",
      [](double ea, double l0, double dx, double dz, double w, double h) {
        const ErnstStiffness e = engine([&] { return ernst_chord_stiffness(ea, l0, dx, dz, w, h); });
        py::dict d;
        d["k"] = e.k;
        d["k_elastic"] = e.k_elastic;
        d["k_sag"] = e.k_sag;
        d["chi"] = e.chi;
        d["note"] = e.note;
        return d;
      },
      py::kw_only(), py::arg("ea"), py::arg("l0"), py::arg("dx"), py::arg("dz"), py::arg("w"), py::arg("h"),
      "Ernst's linear chord stiffness of an elastic catenary: 1/k = l0/EA + (w dx)^2 Lc / (12 T^3), T = h Lc / dx. "
      "dict(k, k_elastic, k_sag, chi, note): chi = w Lc cos^2(alpha) / h the sag parameter (error ~ chi^2 / 5 on "
      "the exact catenary), note = chi > 0.1.");
  py::reinterpret_borrow<py::class_<PyModel>>(mod.attr("Model"))
      .def(
          "deck_camber",
          [](const PyModel& self, Id cs, const std::map<Id, double>& targets, int max_iterations, double tol) {
            ModelPtr model = self.solver_model();
            CamberOptions o;
            o.max_iterations = max_iterations;
            o.tol_rel = tol;
            std::vector<std::pair<Id, double>> t(targets.begin(), targets.end());
            const CamberResults r = engine([&] { return solve_deck_camber(*model, cs, t, o); });
            py::dict camber, dev;
            for (std::size_t i = 0; i < t.size(); ++i) {
              camber[py::int_(t[i].first)] = r.camber[i];
              dev[py::int_(t[i].first)] = r.deviation[i];
            }
            py::dict out;
            out["camber"] = camber;
            out["deviation"] = dev;
            out["iterations"] = r.iterations;
            out["trace"] = r.trace;
            return out;
          },
          py::arg("case"), py::arg("targets"), py::arg("max_iterations") = 30, py::arg("tol") = 1e-9,
          "Fabrication camber: the vertical offsets of the nodes of `targets` {node: target elevation z [m] after "
          "the dead load} such that the deformed profile under case `case` meets them (fixed point "
          "c <- c - (z - target), the structure re-fabricated and re-solved each pass: cable statics with cables, "
          "linear statics otherwise; tol relative to the model size). Returns dict(camber {node: offset, + = "
          "fabricated high}, deviation, iterations, trace).");
  // ---- suspension: form finding of a main cable, assembly, two-stage correction
  const auto form_dict = [](const FormFound& f) {
    const auto v3 = [](const Vec3& v) { return Vec3T(v.x(), v.y(), v.z()); };
    std::vector<Vec3T> nodes, forces;
    for (const Vec3& v : f.nodes) nodes.push_back(v3(v));
    for (const Vec3& v : f.hanger_force) forces.push_back(v3(v));
    py::dict d;
    d["H"] = f.H;
    d["nodes"] = nodes;
    d["segment_l0"] = f.segment_l0;
    d["segment_h"] = f.segment_h;
    d["segment_tmax"] = f.segment_tmax;
    d["hanger_tension"] = f.hanger_tension;
    d["hanger_force"] = forces;
    d["sag"] = f.sag;
    d["saddle_force_left"] = v3(f.saddle_force_left);
    d["saddle_force_right"] = v3(f.saddle_force_right);
    d["total_l0"] = f.total_l0;
    d["iterations"] = f.iterations;
    return d;
  };
  const auto form_spec = [](const Vec3T& left, const Vec3T& right, const std::vector<double>& hanger_x,
                            const std::vector<double>& hanger_load, double ea, double w,
                            std::optional<double> target_sag, std::optional<double> target_h,
                            const std::vector<std::pair<double, double>>& hanger_deck) {
    FormFindSpec f;
    f.left = Vec3(std::get<0>(left), std::get<1>(left), std::get<2>(left));
    f.right = Vec3(std::get<0>(right), std::get<1>(right), std::get<2>(right));
    f.hanger_x = hanger_x;
    f.hanger_load = hanger_load;
    f.EA = ea;
    f.w = w;
    f.target_sag = target_sag;
    f.target_H = target_h;
    for (const auto& [x, z] : hanger_deck) f.hanger_deck.push_back({x, z});
    return f;
  };
  using Anchors = std::vector<std::pair<double, double>>;
  mod.def(
      "form_find_cable",
      [form_spec, form_dict](const Vec3T& left, const Vec3T& right, const std::vector<double>& hanger_x,
                             const std::vector<double>& hanger_load, double ea, double w,
                             std::optional<double> target_sag, std::optional<double> target_h,
                             const Anchors& hanger_deck) {
        const FormFindSpec f = form_spec(left, right, hanger_x, hanger_load, ea, w, target_sag, target_h, hanger_deck);
        return form_dict(engine([&] { return form_find_cable(f); }));
      },
      py::arg("left"), py::arg("right"), py::arg("hanger_x"), py::arg("hanger_load"), py::kw_only(), py::arg("ea"),
      py::arg("w") = 0.0, py::arg("target_sag") = py::none(), py::arg("target_h") = py::none(),
      py::arg("hanger_deck") = Anchors{},
      "Form finding of a suspension main cable between the saddles left and right (x, y, z): hangers attached "
      "at the abscissas hanger_x from the left saddle [m], each carrying the downward load hanger_load [N]; cable "
      "of rigidity ea [N] and weight w [N per unstrained m]; exactly one of target_sag (mid-span sag below the "
      "saddle chord [m]) or target_h (horizontal tension of the first segment [N]). hanger_deck = [(x, z)]: deck "
      "anchor of each hanger in the cable plane, for inclined hangers (default: vertical). Returns dict(H, nodes, "
      "segment_l0, segment_h, segment_tmax, hanger_tension, hanger_force, sag, saddle_force_left, "
      "saddle_force_right, total_l0, iterations): each segment an elastic catenary in equilibrium at its nodes.");
  mod.def(
      "hanger_unstrained_length",
      [](double length, double ea, double tension, double w_axial) {
        return engine([&] { return hanger_unstrained_length(length, ea, tension, w_axial); });
      },
      py::arg("length"), py::arg("ea"), py::arg("tension"), py::arg("w_axial") = 0.0,
      "Unstrained length [m] of a straight hanger of stretched length `length` [m] and axial rigidity ea [N] whose "
      "tension is `tension` [N] at one end and varies by w_axial [N/m] towards the other: length = L0 + (tension "
      "L0 + w_axial L0^2 / 2) / ea. From the foot of a vertical hanger of weight w per unit length w_axial = w, "
      "from its top w_axial = -w; an inclined hanger carries w drop / length along its axis.");

  // The kernel assembles the suspension in a frame model whose cables are bars
  // of kind Cable; the Python model keeps its cables apart: they are taken
  // out of the assembled model and recorded with their unstrained length.
  const auto adopt = [form_dict](PyModel& self, Model& built, const SuspensionModel& sm) {
    std::vector<Id> cable_bars = sm.cable_bars;
    cable_bars.insert(cable_bars.end(), sm.backstay_bars.begin(), sm.backstay_bars.end());
    std::vector<Id> cables;
    for (Id id : cable_bars) {
      const Bar& b = built.bar(id);
      const Material& mat = built.material(b.material);
      const double area = built.section(b.section).props.A;
      CableRec c;
      c.id = self.next_cable++;
      c.n1 = b.ni;
      c.n2 = b.nj;
      c.ea = b.cable && b.cable->ea > 0 ? b.cable->ea : mat.E * area;
      c.w = b.cable && b.cable->w ? *b.cable->w : mat.rho * area * kGravity;
      c.route = 0;
      c.value = b.cable ? b.cable->unstressed_length : 0.0;
      self.cables.push_back(c);
      cables.push_back(c.id);
    }
    std::erase_if(built.bars, [&](const Bar& b) { return b.kind == BarKind::Cable; });
    if (self.active_bars)
      for (Id id : sm.hanger_bars) self.active_bars->insert(id);
    if (self.active_supports)
      for (const Support& s : built.supports) self.active_supports->insert(s.node);
    self.m = built;
    self.m.invalidate_index();
    for (const Node& n : self.m.nodes) self.next_node = std::max(self.next_node, n.id + 1);
    for (const Bar& b : self.m.bars) self.next_bar = std::max(self.next_bar, b.id + 1);
    const std::size_t main = sm.cable_bars.size();
    py::dict d;
    d["form"] = form_dict(sm.form);
    d["cable_nodes"] = sm.cable_nodes;
    d["cables"] = std::vector<Id>(cables.begin(), cables.begin() + static_cast<std::ptrdiff_t>(main));
    d["hanger_bars"] = sm.hanger_bars;
    d["hanger_l0"] = sm.hanger_l0;
    d["backstays"] = std::vector<Id>(cables.begin() + static_cast<std::ptrdiff_t>(main), cables.end());
    d["anchor_nodes"] = sm.anchor_nodes;
    d["deck_nodes"] = sm.deck_nodes;
    d["deck_bars"] = sm.deck_bars;
    return d;
  };
  // The suspension of a dict returned by build_suspension, as the kernel
  // reads it in the solver model (the cables there are bars after the frame bars).
  const auto suspension_of = [](const PyModel& self, const py::dict& d) {
    SuspensionModel sm;
    for (Id c : d["cables"].cast<std::vector<Id>>()) sm.cable_bars.push_back(self.cable_bar_id(c));
    for (Id c : d["backstays"].cast<std::vector<Id>>()) sm.backstay_bars.push_back(self.cable_bar_id(c));
    sm.cable_nodes = d["cable_nodes"].cast<std::vector<Id>>();
    sm.hanger_bars = d["hanger_bars"].cast<std::vector<Id>>();
    sm.deck_nodes = d["deck_nodes"].cast<std::vector<Id>>();
    sm.deck_bars = d["deck_bars"].cast<std::vector<Id>>();
    sm.anchor_nodes = d["anchor_nodes"].cast<std::vector<Id>>();
    return sm;
  };
  const auto response_dict = [](const SuspensionResponse& r) {
    py::dict d;
    d["cable_H"] = r.cable_H;
    d["cable_Tmax"] = r.cable_Tmax;
    d["segment_H"] = r.segment_H;
    d["hanger_tension"] = r.hanger_tension;
    d["deck_my_max"] = r.deck_my_max;
    d["deck_my_min"] = r.deck_my_min;
    d["deck_dz_min"] = r.deck_dz_min;
    d["deck_dz_max"] = r.deck_dz_max;
    d["backstay_H"] = r.backstay_H;
    d["backstay_Tmax"] = r.backstay_Tmax;
    return d;
  };
  const auto backstay_of = [](const py::object& o) -> std::optional<BackstaySpec> {
    if (o.is_none()) return std::nullopt;
    const py::dict d = o.cast<py::dict>();
    BackstaySpec b;
    const Vec3T a = d["anchor"].cast<Vec3T>();
    b.anchor = Vec3(std::get<0>(a), std::get<1>(a), std::get<2>(a));
    if (d.contains("regime")) {
      const std::string r = d["regime"].cast<std::string>();
      if (r != "clamped" && r != "frictionless") throw py::value_error("backstay: regime is 'clamped' or 'frictionless'");
      b.regime = r == "clamped" ? SaddleRegime::Clamped : SaddleRegime::Frictionless;
    }
    if (d.contains("back_h") && !d["back_h"].is_none()) b.back_H = d["back_h"].cast<double>();
    if (d.contains("anchor_node")) b.anchor_node = d["anchor_node"].cast<Id>();
    return b;
  };
  const auto suspension_spec = [form_spec, backstay_of](
                                   const Vec3T& left, const Vec3T& right, const std::vector<double>& hanger_x,
                                   const std::vector<double>& hanger_load, const std::vector<Id>& hanger_deck_nodes,
                                   Id cable_material, Id cable_section, Id hanger_material, Id hanger_section,
                                   std::optional<double> target_sag, std::optional<double> target_h,
                                   const std::string& hanger_kind, Id left_node, Id right_node, bool fix_saddles,
                                   const std::vector<Id>& deck_bars, const py::object& back_left,
                                   const py::object& back_right, const Anchors& hanger_deck,
                                   Id hanger_fit_case, double hanger_weight) {
    if (hanger_kind != "truss" && hanger_kind != "frame")
      throw py::value_error("suspension: hanger_kind is 'truss' or 'frame'");
    SuspensionSpec s;
    s.form = form_spec(left, right, hanger_x, hanger_load, 0.0, 0.0, target_sag, target_h, hanger_deck);
    s.cable_material = cable_material;
    s.cable_section = cable_section;
    s.hanger_material = hanger_material;
    s.hanger_section = hanger_section;
    s.hanger_kind = hanger_kind == "truss" ? BarKind::Truss : BarKind::Frame;
    s.hanger_deck_nodes = hanger_deck_nodes;
    s.left_node = left_node;
    s.right_node = right_node;
    s.fix_saddles = fix_saddles;
    s.deck_bars = deck_bars;
    s.back_left = backstay_of(back_left);
    s.back_right = backstay_of(back_right);
    s.hanger_fit_case = hanger_fit_case;
    s.hanger_weight = hanger_weight;
    return s;
  };
  py::reinterpret_borrow<py::class_<PyModel>>(mod.attr("Model"))
      .def(
          "build_suspension",
          [suspension_spec, adopt](PyModel& self, const Vec3T& left, const Vec3T& right,
                                   const std::vector<double>& hanger_x, const std::vector<double>& hanger_load,
                                   const std::vector<Id>& hanger_deck_nodes, Id cable_material, Id cable_section,
                                   Id hanger_material, Id hanger_section, std::optional<double> target_sag,
                                   std::optional<double> target_h, const std::string& hanger_kind, Id left_node,
                                   Id right_node, bool fix_saddles, const std::vector<Id>& deck_bars,
                                   const py::object& back_left, const py::object& back_right,
                                   const Anchors& hanger_deck, Id hanger_fit_case, double hanger_weight) {
            const SuspensionSpec spec = suspension_spec(
                left, right, hanger_x, hanger_load, hanger_deck_nodes, cable_material, cable_section, hanger_material,
                hanger_section, target_sag, target_h, hanger_kind, left_node, right_node, fix_saddles, deck_bars,
                back_left, back_right, hanger_deck, hanger_fit_case, hanger_weight);
            Model built = self.m;
            const SuspensionModel sm = engine([&] { return build_suspension(built, spec); });
            return adopt(self, built, sm);
          },
          py::arg("left"), py::arg("right"), py::arg("hanger_x"), py::arg("hanger_load"),
          py::arg("hanger_deck_nodes"), py::kw_only(), py::arg("cable_material"), py::arg("cable_section"),
          py::arg("hanger_material"), py::arg("hanger_section"), py::arg("target_sag") = py::none(),
          py::arg("target_h") = py::none(), py::arg("hanger_kind") = "truss", py::arg("left_node") = 0,
          py::arg("right_node") = 0, py::arg("fix_saddles") = true, py::arg("deck_bars") = std::vector<Id>{},
          py::arg("back_left") = py::none(), py::arg("back_right") = py::none(), py::arg("hanger_deck") = Anchors{},
          py::arg("hanger_fit_case") = 0, py::arg("hanger_weight") = 0.0,
          "Form-find a suspension main cable (arguments of keta.form_find_cable; rigidity and weight from "
          "cable_material and cable_section) and assemble it in the model: cable nodes, one cable per segment at "
          "its unstrained length, one hanger bar (hanger_kind 'truss' or 'frame') from each cable node to the deck "
          "node of hanger_deck_nodes. Saddles: existing nodes left_node / right_node, or new nodes fully fixed "
          "(fix_saddles). back_left / back_right = dict(anchor=(x, y, z), regime='clamped'|'frictionless', "
          "back_h=None, anchor_node=0): side spans to an anchorage. hanger_fit_case: the load case that receives "
          "the lack of fit of the hangers, so that each carries its form-found tension at the drawn geometry and "
          "the deck hangs at its drawn level (0: the hangers are bars of their geometric length, and the deck "
          "drops by their stretch); hanger_weight [N/m]: the weight that case applies to a hanger, hanger_load "
          "then being the load at the cable node. Returns dict(form, cable_nodes, cables, hanger_bars, hanger_l0, "
          "backstays, anchor_nodes, deck_nodes, deck_bars): cables and backstays are cable ids of the model, "
          "hanger_l0 the unstrained length of each hanger (keta.hanger_unstrained_length). Solve with "
          "cable_static().")
      .def(
          "solve_uelm",
          [suspension_spec, adopt](PyModel& self, Id dead_case, const Vec3T& left, const Vec3T& right,
                                   const std::vector<double>& hanger_x, const std::vector<double>& hanger_load,
                                   const std::vector<Id>& hanger_deck_nodes, Id cable_material, Id cable_section,
                                   Id hanger_material, Id hanger_section, std::optional<double> target_sag,
                                   std::optional<double> target_h, const std::string& hanger_kind, Id left_node,
                                   Id right_node, bool fix_saddles, const std::vector<Id>& deck_bars, double tol,
                                   int max_iterations, double hanger_weight) {
            if (!self.cables.empty())
              throw py::value_error("solve_uelm: the model already carries cables, which the correction would not see");
            const SuspensionSpec spec = suspension_spec(
                left, right, hanger_x, hanger_load, hanger_deck_nodes, cable_material, cable_section, hanger_material,
                hanger_section, target_sag, target_h, hanger_kind, left_node, right_node, fix_saddles, deck_bars,
                py::none(), py::none(), {}, 0, hanger_weight);
            UelmOptions o;
            o.tuning.tol = tol;
            o.tuning.max_iter = max_iterations;
            Model built = self.m;
            const UelmResult r = engine([&] { return solve_uelm(built, spec, dead_case, o); });
            py::dict d = adopt(self, built, r.sm);
            d["l0"] = r.l0;
            d["deviation_stage1"] = r.deviation_stage1;
            d["deviation"] = r.deviation;
            d["iterations"] = r.iterations;
            return d;
          },
          py::arg("case"), py::arg("left"), py::arg("right"), py::arg("hanger_x"), py::arg("hanger_load"),
          py::arg("hanger_deck_nodes"), py::kw_only(), py::arg("cable_material"), py::arg("cable_section"),
          py::arg("hanger_material"), py::arg("hanger_section"), py::arg("target_sag") = py::none(),
          py::arg("target_h") = py::none(), py::arg("hanger_kind") = "truss", py::arg("left_node") = 0,
          py::arg("right_node") = 0, py::arg("fix_saddles") = true, py::arg("deck_bars") = std::vector<Id>{},
          py::arg("tol") = 1e-9, py::arg("max_iterations") = 50, py::arg("hanger_weight") = 0.0,
          "Two-stage form finding of a suspension bridge (arguments of build_suspension): stage 1 form-finds and "
          "assembles the cable under the hanger loads, the hangers at their unstrained length (their lack of fit "
          "goes to `case`); stage 2 solves the dead case `case` on the assembled "
          "structure — the deck takes part of the load, the cable leaves its profile — and corrects the "
          "unstrained lengths of the main-cable segments until every hanger node is back on the stage-1 profile. "
          "Returns the dict of build_suspension plus l0 (corrected lengths, those of the model's cables), "
          "deviation_stage1 and deviation [m] (largest node deviation before and after), iterations.");
  py::reinterpret_borrow<py::class_<PyModel>>(mod.attr("Model"))
      .def(
          "build_suspension_on_towers",
          [suspension_spec, adopt](PyModel& self, Id dead_case, const Vec3T& left, const Vec3T& right,
                                   const std::vector<double>& hanger_x, const std::vector<double>& hanger_load,
                                   const std::vector<Id>& hanger_deck_nodes, Id left_node, Id right_node,
                                   Id cable_material, Id cable_section, Id hanger_material, Id hanger_section,
                                   std::optional<double> target_sag, std::optional<double> target_h,
                                   const std::string& hanger_kind, const std::vector<Id>& deck_bars,
                                   const py::object& back_left, const py::object& back_right,
                                   const Anchors& hanger_deck, double saddle_tol, int max_iterations,
                                   double hanger_weight) {
            if (!self.cables.empty())
              throw py::value_error(
                  "build_suspension_on_towers: the model already carries cables, which the coupling would not see");
            if (left_node == 0 || right_node == 0)
              throw py::value_error("build_suspension_on_towers: left_node and right_node are the tower-top nodes");
            const SuspensionSpec spec = suspension_spec(
                left, right, hanger_x, hanger_load, hanger_deck_nodes, cable_material, cable_section, hanger_material,
                hanger_section, target_sag, target_h, hanger_kind, left_node, right_node, true, deck_bars, back_left,
                back_right, hanger_deck, 0, hanger_weight);
            SuspensionFlexOptions o;
            o.saddle_tol = saddle_tol;
            o.max_iterations = max_iterations;
            Model built = self.m;
            const SuspensionFlexResult r = engine([&] { return build_suspension_on_towers(built, spec, dead_case, o); });
            py::dict d = adopt(self, built, r.sm);
            d["saddle_disp_left"] = Vec3T(r.saddle_disp_left.x(), r.saddle_disp_left.y(), r.saddle_disp_left.z());
            d["saddle_disp_right"] = Vec3T(r.saddle_disp_right.x(), r.saddle_disp_right.y(), r.saddle_disp_right.z());
            d["saddle_gap"] = r.saddle_gap;
            d["iterations"] = r.iterations;
            return d;
          },
          py::arg("case"), py::arg("left"), py::arg("right"), py::arg("hanger_x"), py::arg("hanger_load"),
          py::arg("hanger_deck_nodes"), py::kw_only(), py::arg("left_node"), py::arg("right_node"),
          py::arg("cable_material"), py::arg("cable_section"), py::arg("hanger_material"),
          py::arg("hanger_section"), py::arg("target_sag") = py::none(), py::arg("target_h") = py::none(),
          py::arg("hanger_kind") = "truss", py::arg("deck_bars") = std::vector<Id>{},
          py::arg("back_left") = py::none(), py::arg("back_right") = py::none(), py::arg("hanger_deck") = Anchors{},
          py::arg("saddle_tol") = 1e-8, py::arg("max_iterations") = 40, py::arg("hanger_weight") = 0.0,
          "Suspension on flexible towers (arguments of build_suspension; left_node and right_node are the existing "
          "tower-top nodes at left and right): the cable is form-found with its saddles where the tower tops stand "
          "under the permanent case `case`, which also receives the lack of fit of the hangers — a fixed point between the form finding and the solve of the whole, "
          "until the saddles move by less than saddle_tol [m]. Returns the dict of build_suspension plus "
          "saddle_disp_left, saddle_disp_right (the tower-top displacements under the permanent case), saddle_gap "
          "(last update) and iterations.")
      .def(
          "suspension_response",
          [suspension_of, response_dict](const PyModel& self, const py::dict& suspension, Id cs) {
            const SuspensionModel sm = suspension_of(self, suspension);
            ModelPtr model = self.solver_model();
            return response_dict(engine([&] { return suspension_response(*model, sm, cs); }));
          },
          py::arg("suspension"), py::arg("case"),
          "Response of a suspension (the dict returned by build_suspension, build_suspension_on_towers or "
          "solve_uelm) under case `case`, by the cable statics of the whole: dict(cable_H, cable_Tmax, segment_H, "
          "hanger_tension, deck_my_max, deck_my_min, deck_dz_min, deck_dz_max, backstay_H, backstay_Tmax) — "
          "horizontal and largest tension of the main cable, tension of each hanger, extreme moment and vertical "
          "displacement of the deck (deck_bars and hanger_deck_nodes of the assembly).")
      .def(
          "suspension_live_response",
          [suspension_of, response_dict](const PyModel& self, const py::dict& suspension, Id dead_case, Id live_case) {
            const SuspensionModel sm = suspension_of(self, suspension);
            ModelPtr model = self.solver_model();
            const SuspensionDelta r = engine([&] { return suspension_live_response(*model, sm, dead_case, live_case); });
            py::dict d;
            d["dead"] = response_dict(r.dead);
            d["total"] = response_dict(r.total);
            d["hanger_delta"] = r.hanger_delta;
            d["cable_H_delta"] = r.cable_H_delta;
            d["deck_my_max_delta"] = r.deck_my_max_delta;
            d["deck_my_min_delta"] = r.deck_my_min_delta;
            return d;
          },
          py::arg("suspension"), py::arg("dead_case"), py::arg("live_case"),
          "Effect of a live load on a suspension: cables do not superpose, so the dead case alone and the dead "
          "and live loads in one case are both solved. Returns dict(dead, total — two dicts of suspension_response "
          "—, hanger_delta, cable_H_delta, deck_my_max_delta, deck_my_min_delta), the differences total - dead.");
  py::reinterpret_borrow<py::class_<PyCableStatic>>(mod.attr("CableStatic"))
      .def(
          "ernst_substitute",
          [](const PyCableStatic& s) {
            const ErnstSubstitute sub = engine([&] { return ernst_substitute(*s.model, *s.r); });
            py::dict axial, chi;
            for (std::size_t i = 0; i < sub.model.bars.size(); ++i) {
              axial[py::int_(sub.model.bars[i].id)] = sub.axial[i];
              if (s.model->bars[i].kind == BarKind::Cable) chi[py::int_(sub.model.bars[i].id)] = sub.cables[i].chi;
            }
            return py::make_tuple(PyModel::from_model(sub.model), axial, chi);
          },
          "The linear substitute of the cable state: (model, axial, chi) — a model whose catenary cables are trusses "
          "of their Ernst stiffness at this state (bars renumbered as in the solver model, cables after the frame "
          "bars), the tension of the state per bar {bar: N} (the prestress of a modal analysis, with it the string "
          "stiffness of the stays) and the sag parameter per cable bar. The trusses lie on the model's geometry.");
}
