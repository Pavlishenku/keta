// keta Python bindings — non-uniform (Vlasov) torsion and flexural-torsional
// buckling of thin-walled bars. Included by keta_py.cpp inside its
// anonymous namespace.
//
//   r = m.warping(case=1, Iw={s: 126e-9}, warping_fixed=[n0])
//   r.displacement(n), r.bar(b)["Mw"]
//   factors = m.lateral_buckling(case=1, Iw={s: 126e-9}).factors

using Vec7 = Eigen::Matrix<double, 7, 1>;

std::array<double, 7> dof7(const Vec7& v) {
  std::array<double, 7> a{};
  for (int i = 0; i < 7; ++i) a[static_cast<std::size_t>(i)] = v(i);
  return a;
}

WarpingSpec warping_spec(const std::map<Id, double>& Iw, const std::vector<Id>& fixed,
                         const std::map<Id, std::pair<double, double>>& offset) {
  WarpingSpec s;
  s.Iw = Iw;
  s.warping_fixed.insert(fixed.begin(), fixed.end());
  for (const auto& [bar, e] : offset) s.load_offset[bar] = Eigen::Vector2d(e.first, e.second);
  return s;
}

struct PyWarping {
  ModelPtr model;
  std::shared_ptr<const WarpingResults> r;
};
struct PyWarpingBuckling {
  ModelPtr model;
  std::shared_ptr<const WarpingBucklingResults> r;
};

void bind_stability(py::module_& mod) {
  auto M = py::reinterpret_borrow<py::class_<PyModel>>(mod.attr("Model"));

  M.def(
      "warping",
      [](PyModel& self, std::optional<Id> c, const std::map<Id, double>& Iw, const std::vector<Id>& fixed,
         const std::map<Id, std::pair<double, double>>& offset, const std::string& theory, int max_iter,
         double tol) {
        ModelPtr model = self.solver_model();
        const Id cid = self.case_or_first(c);
        WarpingOptions o;
        if (theory == "first") o.theory = WarpingTheory::First;
        else if (theory == "second") o.theory = WarpingTheory::Second;
        else throw py::value_error("theory 'first' or 'second'");
        o.max_iterations = max_iter;
        o.tol_rel = tol;
        const WarpingSpec spec = warping_spec(Iw, fixed, offset);
        auto r = engine([&] { return std::make_shared<const WarpingResults>(solve_warping(*model, cid, spec, o)); });
        return PyWarping{model, r};
      },
      py::arg("case") = py::none(), py::arg("Iw") = std::map<Id, double>{},
      py::arg("warping_fixed") = std::vector<Id>{}, py::arg("load_offset") = std::map<Id, std::pair<double, double>>{},
      py::arg("theory") = "first", py::arg("max_iter") = 50, py::arg("tol") = 1e-10,
      "Non-uniform (Vlasov) torsion: a 7th DOF per node, the warping psi = theta_x'. Iw = {section id: C_M "
      "[m^6]} overriding the sections' own Iw (add_section_props(iw=...), the I shapes' computed one); "
      "warping_fixed = nodes whose warping is restrained; load_offset = {bar: (ey, ez)} height of its "
      "span loads above the shear centre. theory 'second' adds the flexural-torsional coupling of the "
      "thin-walled bar (Wagner term N r0^2, moments, load heights) and iterates. Doubly-symmetric sections, "
      "straight bars, no releases / eccentricities / Winkler / CARA VAR (typed refusals).");

  M.def(
      "lateral_buckling",
      [](PyModel& self, std::optional<Id> c, const std::map<Id, double>& Iw, const std::vector<Id>& fixed,
         const std::map<Id, std::pair<double, double>>& offset, int n_modes) {
        ModelPtr model = self.solver_model();
        const Id cid = self.case_or_first(c);
        const WarpingSpec spec = warping_spec(Iw, fixed, offset);
        auto r = engine([&] {
          return std::make_shared<const WarpingBucklingResults>(warping_buckling(*model, cid, spec, n_modes));
        });
        return PyWarpingBuckling{model, r};
      },
      py::arg("case") = py::none(), py::arg("Iw") = std::map<Id, double>{},
      py::arg("warping_fixed") = std::vector<Id>{}, py::arg("load_offset") = std::map<Id, std::pair<double, double>>{},
      py::arg("n_modes") = 3,
      "Linear flexural-torsional (lateral-torsional) buckling about the first-order state of the case: "
      "(K + lambda Kg) x = 0 with Kg from N, My, Mz and the load heights; factors scale that case. The "
      "coupling is a discretization (Hermite cubics): mesh the members (h^2 to h^4 convergence).");

  M.def(
      "buckling_cte_var",
      [](PyModel& self, Id cte, Id var, int n_modes, bool tension_stiffening, bool bisection, bool exact_stability) {
        ModelPtr model = self.solver_model();
        py::dict d;
        if (bisection) {
          BucklingCteVarBisectOptions o;
          o.exact_stability = exact_stability;
          o.tension_stiffening = tension_stiffening;
          d["factors"] = std::vector<double>{engine([&] { return buckling_cte_var_bisection(*model, cte, var, o); })};
          return d;
        }
        BucklingCteVarOptions o;
        o.tension_stiffening = tension_stiffening;
        const BucklingCteVar r = engine([&] { return solve_buckling_cte_var(*model, cte, var, n_modes, o); });
        d["factors"] = r.load_factors;
        d["residuals"] = r.residuals;
        std::vector<Id> ids;
        for (const Node& n : model->nodes) ids.push_back(n.id);
        d["nodes"] = ids;
        std::vector<std::vector<Dof6>> modes;
        for (const auto& md : r.modes) {
          std::vector<Dof6> v;
          for (const Vec6& x : md) v.push_back(dof6(x));
          modes.push_back(v);
        }
        d["modes"] = modes;
        d["N_cte"] = r.N_cte;
        d["N_var"] = r.N_var;
        d["gen_stiffness"] = r.gen_stiffness;
        d["participation"] = r.participation;
        d["total_kg"] = r.total_kg;
        return d;
      },
      py::arg("cte"), py::arg("var"), py::arg("n_modes") = 2, py::arg("tension_stiffening") = false,
      py::arg("bisection") = false, py::arg("exact_stability") = true,
      "Buckling on a fixed state plus an amplified one: K + Kg(N_cte) + lambda Kg(N_var), the axial fields "
      "frozen from load cases cte (0 = none) and var. Returns dict(factors, residuals, nodes, modes (unit "
      "crest, per node in model order), N_cte, N_var, gen_stiffness, participation, total_kg). bisection=True "
      "gives the critical factor alone by LDLT inertia bisection, with exact_stability member tangents (one "
      "element per member then suffices). Tension does not stiffen unless tension_stiffening.");

  py::class_<PyWarping>(mod, "WarpingResult", "Static state with non-uniform torsion (7 DOFs per node).")
      .def("displacement",
           [](const PyWarping& s, Id node) { return dof7(s.r->node_disp.at(node_pos(*s.model, node))); },
           "(ux, uy, uz, rx, ry, rz, psi) at a node, psi = theta_x' [1/m].")
      .def("bar",
           [](const PyWarping& s, Id bar) {
             for (const WarpingBarResult& b : s.r->bars)
               if (b.bar == bar) {
                 py::dict d;
                 d["N"] = b.N;
                 d["MT"] = b.MT;
                 d["MT1"] = b.MT1;
                 d["MT2"] = b.MT2;
                 d["MT3"] = b.MT3;
                 d["Mw"] = b.Mw;
                 d["end_forces"] = std::vector<double>(b.end_forces.data(), b.end_forces.data() + 14);
                 return d;
               }
             throw py::key_error("no bar " + std::to_string(bar));
           },
           "dict(N, MT, MT1 (Saint-Venant), MT2 (warping), MT3 (Wagner), Mw (bimoment) at ends (i, j), "
           "end_forces = 12 local member end forces + the two bimoments).")
      .def("reactions",
           [](const PyWarping& s) {
             std::vector<std::pair<Id, std::array<double, 7>>> out;
             for (const auto& [n, v] : s.r->reactions) out.emplace_back(n, dof7(v));
             return out;
           })
      .def_property_readonly("iterations", [](const PyWarping& s) { return s.r->iterations; });

  py::class_<PyWarpingBuckling>(mod, "WarpingBuckling", "Flexural-torsional buckling factors and modes.")
      .def_property_readonly("factors", [](const PyWarpingBuckling& s) { return s.r->load_factors; })
      .def(
          "mode",
          [](const PyWarpingBuckling& s, int i, Id node) {
            if (i < 0 || static_cast<std::size_t>(i) >= s.r->modes.size()) throw py::index_error("mode index");
            return dof7(s.r->modes[static_cast<std::size_t>(i)].at(node_pos(*s.model, node)));
          },
          py::arg("i"), py::arg("node"), "Unit-crest mode i at a node (ux, uy, uz, rx, ry, rz, psi).");
}
