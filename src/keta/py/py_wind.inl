// keta Python bindings — wind: profile, turbulence spectra,
// admittance, spectral and temporal buffeting, flutter, equivalent static
// loads, concomitance polyhedron, the second-order wind loop, vortex
// shedding. Included by keta_py.cpp inside its anonymous namespace.
//
//   spec = dict(vb=26, z0=0.05, B=12, aero=dict(alpha=[-0.2, 0.2], cd=[0.6, 0.6],
//               cl=[-0.8, 0.8], cm=[-0.1, 0.1]), zeta=0.01, n_modes=6)
//   r = m.buffeting(spec, [(node, 2)])     # [dict(sigma, sigma_b, sigma_r, nu, g)]

using Matrix = std::vector<std::vector<double>>;
Eigen::MatrixXd to_eigen(const Matrix& a) {
  const Eigen::Index n = static_cast<Eigen::Index>(a.size());
  Eigen::MatrixXd m(n, n == 0 ? 0 : static_cast<Eigen::Index>(a.front().size()));
  for (Eigen::Index i = 0; i < n; ++i) {
    if (a[static_cast<std::size_t>(i)].size() != static_cast<std::size_t>(m.cols())) throw py::value_error("matrix: ragged rows");
    for (Eigen::Index j = 0; j < m.cols(); ++j) m(i, j) = a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
  }
  return m;
}
Matrix from_eigen(const Eigen::MatrixXd& m) {
  Matrix a(static_cast<std::size_t>(m.rows()), std::vector<double>(static_cast<std::size_t>(m.cols())));
  for (Eigen::Index i = 0; i < m.rows(); ++i)
    for (Eigen::Index j = 0; j < m.cols(); ++j) a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = m(i, j);
  return a;
}

Admittance admittance_of(const py::handle& h) {
  Admittance a;
  if (py::isinstance<py::str>(h)) {
    const std::string s = h.cast<std::string>();
    if (s == "quasi_steady") a.kind = Admittance::Kind::QuasiSteady;
    else if (s == "sears") a.kind = Admittance::Kind::Sears;
    else if (s == "liepmann") a.kind = Admittance::Kind::Liepmann;
    else throw py::value_error("admittance 'quasi_steady', 'sears', 'liepmann' or a table [(K, chi2)]");
  } else {
    a.kind = Admittance::Kind::Table;
    a.table = h.cast<std::vector<std::pair<double, double>>>();
  }
  return a;
}

BuffetingSpec spec_of(const py::dict& d) {
  static const std::set<std::string> known = {"vb", "z0", "zmin", "rho", "co", "kI", "spectrum", "rv", "rw", "lv",
                                              "lw", "coherence", "wind_dir", "B", "aero", "alpha0", "admittance",
                                              "aeroelastic", "zeta", "n_log", "per_mode", "duration", "n_modes"};
  for (const auto& [k, v] : d)
    if (!known.count(py::str(k))) throw py::value_error("wind spec: unknown key '" + std::string(py::str(k)) + "'");
  BuffetingSpec s;
  WindProfile& p = s.turb.profile;
  p.vb = opt_get<double>(d, "vb", p.vb);
  p.z0 = opt_get<double>(d, "z0", p.z0);
  p.zmin = opt_get<double>(d, "zmin", p.zmin);
  p.rho = opt_get<double>(d, "rho", p.rho);
  p.co = opt_get<double>(d, "co", p.co);
  p.kI = opt_get<double>(d, "kI", p.kI);
  const std::string form = opt_get<std::string>(d, "spectrum", "EN1991");
  if (form == "EN1991") s.turb.form = SpectrumForm::EN1991;
  else if (form == "von_karman") s.turb.form = SpectrumForm::VonKarman;
  else throw py::value_error("spectrum 'EN1991' or 'von_karman'");
  s.turb.rv = opt_get<double>(d, "rv", s.turb.rv);
  s.turb.rw = opt_get<double>(d, "rw", s.turb.rw);
  s.turb.lv = opt_get<double>(d, "lv", s.turb.lv);
  s.turb.lw = opt_get<double>(d, "lw", s.turb.lw);
  if (d.contains("coherence")) {
    const auto c = d["coherence"].cast<std::array<std::array<double, 2>, 3>>();
    s.turb.C = c;
  }
  const auto wd = opt_get<std::array<double, 3>>(d, "wind_dir", {1, 0, 0});
  s.wind_dir = Vec3(wd[0], wd[1], wd[2]);
  s.B = opt_get<double>(d, "B", 1.0);
  if (d.contains("aero")) {
    const py::dict a = d["aero"].cast<py::dict>();
    if (a.contains("alpha")) {
      s.aero.alpha = a["alpha"].cast<std::vector<double>>();
      s.aero.cd = a["cd"].cast<std::vector<double>>();
      s.aero.cl = a["cl"].cast<std::vector<double>>();
      s.aero.cm = a["cm"].cast<std::vector<double>>();
    } else {
      AeroSection sec;
      sec.cd = opt_get<double>(a, "cd", 0.0);
      sec.cl = opt_get<double>(a, "cl", 0.0);
      sec.cm = opt_get<double>(a, "cm", 0.0);
      sec.dcl = opt_get<double>(a, "dcl", 0.0);
      sec.dcm = opt_get<double>(a, "dcm", 0.0);
      s.aero = AeroTable::from_section(sec);
    }
  } else {
    s.aero = AeroTable::from_section(AeroSection{});
  }
  s.alpha0 = opt_get<double>(d, "alpha0", 0.0);
  if (d.contains("admittance")) {
    const py::handle h = d["admittance"];
    if (py::isinstance<py::list>(h) && py::len(h) == 3 && !py::isinstance<py::tuple>(py::list(h.cast<py::list>())[0])) {
      std::size_t i = 0;
      for (const py::handle x : h) s.chi[i++] = admittance_of(x);
    } else {
      const Admittance a = admittance_of(h);
      s.chi = {a, a, a};
    }
  }
  s.aeroelastic = opt_get<bool>(d, "aeroelastic", false);
  s.damping.zeta = opt_get<double>(d, "zeta", 0.01);
  s.n_log = opt_get<int>(d, "n_log", s.n_log);
  s.per_mode = opt_get<int>(d, "per_mode", s.per_mode);
  s.duration = opt_get<double>(d, "duration", s.duration);
  return s;
}
std::vector<BuffetingRequest> requests_of(const std::vector<std::pair<Id, int>>& r) {
  std::vector<BuffetingRequest> out;
  for (const auto& [n, c] : r) out.push_back({n, c});
  return out;
}
py::dict stats_dict(const ResponseStats& s) {
  py::dict d;
  d["sigma"] = s.sigma;
  d["sigma_b"] = s.sigma_b;
  d["sigma_r"] = s.sigma_r;
  d["nu"] = s.nu;
  d["g"] = s.peak_factor;
  return d;
}

void bind_wind(py::module_& mod) {
  auto cls = py::reinterpret_borrow<py::class_<PyModel>>(mod.attr("Model"));

  py::class_<WindProfile>(mod, "WindProfile", "EN 1991-1-4 mean wind profile (log law).")
      .def(py::init([](double vb, double z0, double zmin, double rho, double co, double kI) {
             WindProfile p;
             p.vb = vb;
             p.z0 = z0;
             p.zmin = zmin;
             p.rho = rho;
             p.co = co;
             p.kI = kI;
             return p;
           }),
           py::arg("vb") = 26.0, py::arg("z0") = 0.05, py::arg("zmin") = 2.0, py::arg("rho") = 1.225, py::arg("co") = 1.0,
           py::arg("kI") = 1.0)
      .def_readwrite("vb", &WindProfile::vb)
      .def_readwrite("z0", &WindProfile::z0)
      .def_readwrite("zmin", &WindProfile::zmin)
      .def_readwrite("rho", &WindProfile::rho)
      .def("kr", &WindProfile::kr, "Terrain factor 0.19 (z0 / 0.05)^0.07.")
      .def("vm", &WindProfile::vm, py::arg("z"), "Mean speed at height z [m/s].")
      .def("iv", &WindProfile::iv, py::arg("z"), "Turbulence intensity at z.")
      .def("qp", &WindProfile::qp, py::arg("z"), "Peak velocity pressure at z [Pa].");

  mod.def("wind_spectrum",
          [](const std::string& form, int comp, double n, double U, double sigma, double L) {
            const SpectrumForm f = form == "EN1991" ? SpectrumForm::EN1991
                                   : form == "von_karman" ? SpectrumForm::VonKarman
                                                          : throw py::value_error("spectrum 'EN1991' or 'von_karman'");
            return engine([&] { return auto_spectrum(f, comp, n, U, sigma, L); });
          },
          py::arg("form"), py::arg("component"), py::arg("n"), py::arg("U"), py::arg("sigma"), py::arg("L"),
          "One-sided spectrum S(n) of component 0 (u), 1 (v), 2 (w) [(m/s)^2/Hz]; integrates to sigma^2.");
  mod.def("davenport_coherence", &davenport_coherence, py::arg("n"), py::arg("dy"), py::arg("dz"), py::arg("Uj"),
          py::arg("Uk"), py::arg("Cy"), py::arg("Cz"));
  mod.def("theodorsen", [](double K) { return engine([&] { return theodorsen_exact(K); }); }, py::arg("K"),
          "Theodorsen's C(K) (exact, Hankel functions), K = omega b / U.");
  mod.def("sears", [](double K) { return engine([&] { return sears_function(K); }); }, py::arg("K"));
  mod.def("admittance",
          [](const py::object& kind, double K) { return admittance_of(kind).chi2(K); }, py::arg("kind"), py::arg("K"),
          "chi^2(K) of 'quasi_steady', 'sears', 'liepmann' or a table [(K, chi2)].");
  mod.def("vortex_critical_speed", &vortex_critical_speed, py::arg("natural_freq"), py::arg("width"),
          py::arg("strouhal") = 0.18);
  mod.def("scruton_number", &scruton_number, py::arg("mass_per_length"), py::arg("damping_ratio"), py::arg("rho"),
          py::arg("width"));
  mod.def("vortex_amplitude_ratio", &vortex_amplitude_ratio, py::arg("strouhal"), py::arg("scruton"), py::arg("K") = 0.13,
          py::arg("Kw") = 0.6, py::arg("clat") = 0.7);

  cls.def(
      "add_mean_wind_case",
      [](PyModel& self, const py::dict& spec, Id cs) {
        const BuffetingSpec s = spec_of(spec);
        const LoadCase lc = engine([&] { return mean_buffeting_case(*self.solver_model(), s, cs, "mean wind"); });
        for (const Load& l : lc.loads) self.case_mut(cs).loads.push_back(l);
      },
      py::arg("spec"), py::arg("case") = 1,
      "Add the mean wind of `spec` (pressure rho U(z)^2 / 2, drag, lift, moment of the aerodynamic table at the mean "
      "incidence) to a load case.");
  cls.def(
      "buffeting",
      [](const PyModel& self, const py::dict& spec, const std::vector<std::pair<Id, int>>& requests, int n_modes) {
        ModelPtr model = self.solver_model();
        const BuffetingSpec s = spec_of(spec);
        ModalOptions mo;
        mo.n_modes = n_modes;
        const MultiBuffetingResults r = engine([&] {
          py::gil_scoped_release nogil;
          return spectral_buffeting_multi(*model, solve_modal(*model, mo), s, requests_of(requests));
        });
        py::list out;
        for (const ResponseStats& st : r.stats) out.append(stats_dict(st));
        py::dict res;
        res["responses"] = out;
        res["modal_cov"] = from_eigen(r.modal_cov);
        return res;
      },
      py::arg("spec"), py::arg("requests"), py::arg("n_modes") = 6,
      "Spectral buffeting (u and w, drag / lift / moment, coherence, admittance, exact cross-modal terms; "
      "aeroelastic damping and stiffness when spec['aeroelastic']). requests: [(node, comp)]. Returns "
      "dict(responses=[dict(sigma, sigma_b, sigma_r, nu, g)], modal_cov). Spec keys: vb, z0, zmin, rho, co, kI, "
      "spectrum ('EN1991' | 'von_karman'), rv, rw, lv, lw, coherence [[Cy, Cz] x 3], wind_dir, B, aero "
      "(dict(alpha, cd, cl, cm) table or dict(cd, cl, cm, dcl, dcm)), alpha0, admittance, aeroelastic, zeta, n_log, "
      "per_mode, duration.");
  cls.def(
      "temporal_buffeting",
      [](const PyModel& self, const py::dict& spec, const std::vector<std::pair<Id, int>>& requests, int n_modes,
         double n_max, int n_freq, std::uint64_t seed, int periods) {
        ModelPtr model = self.solver_model();
        const BuffetingSpec s = spec_of(spec);
        ModalOptions mo;
        mo.n_modes = n_modes;
        const TemporalBuffeting t = engine([&] {
          py::gil_scoped_release nogil;
          return temporal_buffeting(*model, solve_modal(*model, mo), s, requests_of(requests), n_max, n_freq, seed, periods);
        });
        py::list out;
        for (const TemporalResponse& r : t.responses) {
          py::dict d;
          d["sigma"] = r.sigma;
          d["max"] = r.max;
          d["min"] = r.min;
          d["series"] = r.series;
          out.append(d);
        }
        return py::make_tuple(t.dt, out);
      },
      py::arg("spec"), py::arg("requests"), py::arg("n_modes") = 6, py::arg("n_max") = 5.0, py::arg("n_freq") = 1024,
      py::arg("seed") = 1, py::arg("periods") = 2,
      "Time-domain buffeting: wind synthesized at the strips (seeded), quasi-steady forces, Newmark on the modes. "
      "Returns (dt, [dict(sigma, max, min, series)]).");
  cls.def(
      "flutter",
      [](const PyModel& self, double b, double a, int n_modes, std::array<double, 3> wind_dir, double u_min, double u_max,
         int n_sweep, double zeta, std::optional<std::vector<std::array<double, 9>>> scanlan, double rho) {
        ModelPtr model = self.solver_model();
        ModalOptions mo;
        mo.n_modes = n_modes;
        const ModalResults modal = engine([&] { return solve_modal(*model, mo); });
        FlutterStrip proto;
        if (scanlan) {
          proto.kind = FlutterStrip::Kind::Scanlan;
          proto.B = 2 * b;
          proto.der.rows = *scanlan;
        } else {
          proto.kind = FlutterStrip::Kind::Theodorsen;
          proto.b = b;
          proto.a = a;
        }
        MultimodalFlutter f;
        f.omega = modal.omegas;
        f.zeta.assign(modal.omegas.size(), zeta);
        f.rho = rho;
        f.strips = engine([&] { return flutter_strips(*model, modal, Vec3(wind_dir[0], wind_dir[1], wind_dir[2]), proto); });
        const FlutterResult r = engine([&] { return multimodal_flutter(f, u_min, u_max, n_sweep); });
        return py::make_tuple(r.found, r.U_flutter, r.omega_flutter);
      },
      py::arg("half_chord"), py::arg("a") = 0.0, py::arg("n_modes") = 6,
      py::arg("wind_dir") = std::array<double, 3>{1, 0, 0}, py::arg("u_min") = 1.0, py::arg("u_max") = 150.0,
      py::arg("n_sweep") = 80, py::arg("zeta") = 0.0, py::arg("scanlan") = py::none(), py::arg("rho") = 1.225,
      "Multimodal flutter on the model's modes: Theodorsen flat plate (half chord b, axis a half-chords aft of "
      "mid-chord) or Scanlan derivatives rows [(K, H1..H4, A1..A4)] (K = B omega / U, B = 2 b). Returns (found, "
      "U_flutter, omega_flutter).");
  cls.def(
      "wind_loop",
      [](const PyModel& self, const py::dict& spec, std::pair<Id, int> target, std::optional<Id> permanent, int n_modes,
         int max_iterations, double tol) {
        ModelPtr model = self.solver_model();
        const BuffetingSpec s = spec_of(spec);
        const WindLoopResult r = engine([&] {
          return wind_nonlinear_loop(*model, permanent.value_or(0), s, {target.first, target.second}, n_modes,
                                     max_iterations, tol);
        });
        py::dict d;
        d["sigma"] = r.sigma;
        d["g"] = r.peak_factor;
        d["extreme"] = r.extreme;
        d["iterations"] = r.iterations;
        d["converged"] = r.converged;
        d["frequencies_hz"] = r.frequencies_hz;
        return d;
      },
      py::arg("spec"), py::arg("target"), py::arg("permanent_case") = py::none(), py::arg("n_modes") = 6,
      py::arg("max_iterations") = 20, py::arg("tol") = 1e-6,
      "Second-order wind loop: static (permanent + mean + LRC peak) -> axial forces -> modes about them -> "
      "buffeting, to convergence of the target's sigma and extreme.");
  cls.def(
      "lrc_peak_case",
      [](PyModel& self, const py::dict& spec, std::pair<Id, int> target, int n_modes, Id cs) {
        ModelPtr model = self.solver_model();
        const BuffetingSpec s = spec_of(spec);
        ModalOptions mo;
        mo.n_modes = n_modes;
        const ModalResults modal = engine([&] { return solve_modal(*model, mo); });
        const MultiBuffetingResults r =
            engine([&] { return spectral_buffeting_multi(*model, modal, s, {{target.first, target.second}}); });
        const Assembled asm_ = assemble_elements(*model);
        const int f = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(model->node_index(target.first), target.second))];
        if (f < 0) throw py::value_error("lrc_peak_case: the target degree of freedom is restrained");
        Eigen::VectorXd psi(static_cast<Eigen::Index>(modal.omegas.size()));
        for (std::size_t i = 0; i < modal.omegas.size(); ++i) psi(static_cast<Eigen::Index>(i)) = modal.modes_free[i](f);
        const Eigen::VectorXd qe = lrc_modal(r.modal_cov, psi, r.stats[0].peak_factor);
        const LoadCase lc = modal_inertia_case(*model, modal, qe, cs, "LRC peak");
        for (const Load& l : lc.loads) self.case_mut(cs).loads.push_back(l);
        return py::make_tuple(r.stats[0].sigma, r.stats[0].peak_factor);
      },
      py::arg("spec"), py::arg("target"), py::arg("n_modes") = 6, py::arg("case") = 1,
      "Add to a case the LRC equivalent static load of the buffeting peak of `target` (modal inertia loads: a "
      "linear static solve gives g sigma there). Returns (sigma, g).");
  mod.def("lrc_loads",
          [](const Matrix& C, const std::vector<double>& I, double g) {
            const Eigen::VectorXd v = engine([&] {
              return lrc_loads(to_eigen(C), Eigen::Map<const Eigen::VectorXd>(I.data(), static_cast<Eigen::Index>(I.size())), g);
            });
            return std::vector<double>(v.data(), v.data() + v.size());
          },
          py::arg("cov"), py::arg("influence"), py::arg("g"), "Kasperski's LRC load g C I / sqrt(I^T C I).");
  mod.def("equiprobable_polyhedron",
          [](const std::vector<double>& m, const Matrix& S, double g) {
            const Eigen::VectorXd mm = Eigen::Map<const Eigen::VectorXd>(m.data(), static_cast<Eigen::Index>(m.size()));
            std::vector<std::vector<double>> out;
            for (const Eigen::VectorXd& v : engine([&] { return equiprobable_polyhedron(mm, to_eigen(S), g); }))
              out.emplace_back(v.data(), v.data() + v.size());
            return out;
          },
          py::arg("mean"), py::arg("cov"), py::arg("g"),
          "Capra-Davidovici vertices circumscribing the equiprobability ellipse (2 effects: 8) / ellipsoid (3: 24).");
  mod.def("scanlan_from_theodorsen",
          [](double a, const std::vector<double>& Ks) { return engine([&] { return scanlan_from_theodorsen(a, Ks); }).rows; },
          py::arg("a"), py::arg("K"), "Flat-plate Scanlan derivatives rows (K, H1..H4, A1..A4).");
}
