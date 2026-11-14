// keta Python bindings — time-dependent laws: EN 1992-1-1 creep
// and shrinkage, the coded concrete laws (BPEL, BHP, CEB, EC2, fib, B3,
// IP2 ...), relaxation (EC2, BPEL, IP2), deferred losses, the hereditary
// (Volterra) integral and Kelvin-chain fits. Included by keta_py.cpp inside
// its anonymous namespace.
//
//   law = keta.concrete_law("EC2", fc=35e6, rh=70, h_mm=300, cement="N")
//   law.phi(10000, 28), law.eps_sh(10000, 3), law.E(28)
//   eps = keta.volterra_strain(law.compliance, t, sigma)

ConcreteTime concrete_time(double fck, double rh, double h0, const std::string& cement, double t_shrink0,
                           bool autogenous) {
  ConcreteTime c;
  c.fck = fck;
  c.rh = rh;
  c.h0 = h0;
  if (cement == "S") c.cement_s = 0.38;
  else if (cement == "N") c.cement_s = 0.25;
  else if (cement == "R") c.cement_s = 0.20;
  else throw py::value_error("cement class 'S', 'N' or 'R'");
  c.t_shrink0 = t_shrink0;
  c.autogenous = autogenous;
  return c;
}

struct PyMatLaw {
  MatLaw law;
};

ComplianceFn compliance_of(const py::object& J) {
  if (py::isinstance<PyMatLaw>(J)) return matlaw_compliance(J.cast<const PyMatLaw&>().law);
  if (!PyCallable_Check(J.ptr())) throw py::type_error("J: a concrete_law or a callable J(t, t0)");
  return [J](double t, double t0) {
    py::gil_scoped_acquire g;
    return J(t, t0).cast<double>();
  };
}

void bind_creep(py::module_& mod) {
  py::class_<PyMatLaw>(mod, "ConcreteLaw", "A coded concrete law: modulus, strength, creep and shrinkage vs age.")
      .def("E", [](const PyMatLaw& s, double t) { return engine([&] { return s.law.E(t); }); }, py::arg("t"),
           "Young's modulus at age t [days].")
      .def("Ev", [](const PyMatLaw& s, double t) { return engine([&] { return s.law.Ev(t); }); }, py::arg("t"),
           "BPEL / BHP deferred modulus 3700 fcj^(1/3) MPa.")
      .def("E_ref", [](const PyMatLaw& s, double t0) { return engine([&] { return s.law.E_ref(t0); }); },
           py::arg("t0") = 28.0, "The modulus phi is referenced to.")
      .def("fc_mean", [](const PyMatLaw& s, double t) { return engine([&] { return s.law.fc_mean(t); }); },
           py::arg("t"))
      .def("phi", [](const PyMatLaw& s, double t, double t0) { return engine([&] { return s.law.phi(t, t0); }); },
           py::arg("t"), py::arg("t0"), "Creep coefficient phi(t, t0), ages in days.")
      .def("eps_sh", [](const PyMatLaw& s, double t, double ts) { return engine([&] { return s.law.eps_sh(t, ts); }); },
           py::arg("t"), py::arg("ts"), "Shrinkage at age t, drying from ts (negative = shortening).")
      .def("n_equiv", [](const PyMatLaw& s) { return s.law.n_equiv(); })
      .def(
          "compliance",
          [](const PyMatLaw& s, double t, double t0) { return engine([&] { return matlaw_compliance(s.law)(t, t0); }); },
          py::arg("t"), py::arg("t0"), "J(t, t0) = 1/E(t0) + phi(t, t0)/E_ref(t0) [1/Pa].");

  mod.def(
      "concrete_law", [](const std::string& code, const py::kwargs& kw) { return PyMatLaw{make_law(code, kw)}; },
      py::arg("code"),
      "A coded concrete law: code DIRECT, BPEL, BHP1, BHP2, BPEL_BL, CEB, EC2, EC2_BHP1, EC2_BHP2, FIB2010, "
      "ACI_B3, EC2_2027, IP2, MC78 or COURBON; parameters fc [Pa], rh [%], h_mm, rm_cm, cement, ros, ts_dry, "
      "es, er, theta, kb, k1_consistency, phi_f2, silica, fib_secant, e_direct (as Model.set_concrete_law).");

  mod.def(
      "ec2_phi",
      [](double t, double t0, double fck, double rh, double h0, const std::string& cement) {
        return engine([&] { return ec2_phi(t, t0, concrete_time(fck, rh, h0, cement, 3.0, true)); });
      },
      py::arg("t"), py::arg("t0"), py::arg("fck"), py::arg("rh") = 70.0, py::arg("h0") = 0.3, py::arg("cement") = "N",
      "EN 1992-1-1 annex B creep coefficient phi(t, t0); fck [Pa], h0 = 2 Ac / u [m], ages in days.");
  mod.def(
      "ec2_shrinkage",
      [](double t, double fck, double rh, double h0, const std::string& cement, double ts, bool autogenous) {
        return engine([&] { return ec2_shrinkage(t, concrete_time(fck, rh, h0, cement, ts, autogenous)); });
      },
      py::arg("t"), py::arg("fck"), py::arg("rh") = 70.0, py::arg("h0") = 0.3, py::arg("cement") = "N",
      py::arg("ts") = 3.0, py::arg("autogenous") = true,
      "EN 1992-1-1 §3.1.4 total shrinkage (drying from ts + autogenous) at age t; negative = shortening.");
  mod.def("ec2_ecm", &ec2_ecm, py::arg("fck"), "Ecm = 22000 (fcm / 10)^0.3 MPa, fcm = fck + 8 MPa [Pa].");
  mod.def("aaem_modulus", &aaem_modulus, py::arg("E"), py::arg("dphi"), py::arg("chi") = 0.8,
          "Age-adjusted effective modulus E / (1 + chi dphi).");
  mod.def("nonlinear_creep_factor", &nonlinear_creep_factor, py::arg("sigma_c"), py::arg("fcm_t0"),
          "EN 1992-1-1 §3.1.4(4): exp(1.5 (k_sigma - 0.45)), k_sigma = max(0.45, |sigma_c| / fcm(t0)).");
  mod.def(
      "ec2_relaxation_ratio",
      [](double t_hours, double sigma_pi, double fpk, int cls, double rho1000, std::optional<double> seuil,
         bool capped) {
        return engine([&] {
          return capped ? relaxation_ratio_seuil(t_hours, sigma_pi, fpk, cls, rho1000, seuil)
                        : ec2_relaxation_ratio(t_hours, sigma_pi, fpk, cls, rho1000);
        });
      },
      py::arg("t_hours"), py::arg("sigma_pi"), py::arg("fpk"), py::arg("cls") = 2, py::arg("rho1000") = 2.5,
      py::arg("seuil") = py::none(), py::arg("capped") = false,
      "EN 1992-1-1 eq. 3.28-3.30 relaxation loss ratio d_sigma_pr / sigma_pi at t hours (class 1, 2 or 3); "
      "capped=True applies the time threshold `seuil` [h] (0 = 5e5 h).");
  mod.def("bpel_relaxation", &bpel_relaxation, py::arg("t_days"), py::arg("sigma"), py::arg("fprg"),
          py::arg("rho1000"), "BPEL 91 annex 2 pure relaxation loss after t days (unit of sigma).");
  mod.def("ip2_relaxation", &ip2_relaxation, py::arg("t_hours"), py::arg("sigma_pi"), py::arg("rg"),
          py::arg("rho1000"), py::arg("rho3000"), "IP2 art. 9.2 pure relaxation loss (unit of sigma_pi).");
  mod.def("ec2_deferred_loss", &ec2_deferred_loss, py::arg("eps_cs"), py::arg("Ep"), py::arg("d_sigma_pr"),
          py::arg("Ecm"), py::arg("phi"), py::arg("sigma_c"), py::arg("Ap"), py::arg("Ac"), py::arg("Ic"),
          py::arg("zcp"),
          "EN 1992-1-1 eq. 5.46 deferred loss (|eps_cs| Ep + 0.8 d_sigma_pr + Ep/Ecm phi sigma_c) / (1 + Ep Ap / "
          "(Ecm Ac) (1 + Ac zcp^2 / Ic)(1 + 0.8 phi)); sigma_c compression-positive at the tendon.");

  // ---- lump-sum deferred losses (keta/prestress_tools.h) ----
  const auto lump_dict = [](const LumpSumLoss& l) {
    py::dict d;
    d["elastic"] = l.elastic;
    d["shrinkage"] = l.shrinkage;
    d["creep"] = l.creep;
    d["relaxation"] = l.relaxation;
    d["viscoelastic"] = l.viscoelastic;
    d["total"] = l.total();
    return d;
  };
  const auto lump_options = [](double n, const std::optional<double>& eps_r, const std::optional<double>& mu0) {
    LumpSumOptions o;
    o.n = n;
    if (eps_r) {
      if (!(*eps_r >= 0)) throw py::value_error("lump sum: eps_r >= 0 (a shrinkage strain, as a positive value)");
      o.eps_r = *eps_r;
    }
    if (mu0) {
      if (!(*mu0 >= 0)) throw py::value_error("lump sum: mu0 >= 0");
      o.mu0 = *mu0;
    }
    return o;
  };
  mod.def(
      "ip2_lump_sum",
      [lump_dict, lump_options](double sigma_b, double sigma_pi, double fprg, double Ep, double rh,
                                double final_relaxation, double n, std::optional<double> eps_r) {
        const LumpSumOptions o = lump_options(n, eps_r, std::nullopt);
        return lump_dict(engine([&] { return ip2_lump_sum(sigma_b, sigma_pi, fprg, Ep, rh, final_relaxation, o); }));
      },
      py::arg("sigma_b"), py::arg("sigma_pi"), py::arg("fprg"), py::arg("Ep"), py::arg("rh"),
      py::arg("final_relaxation"), py::arg("n") = 6.0, py::arg("eps_r") = py::none(),
      "Lump-sum deferred losses of the IP2 rules from the final concrete stress sigma_b at the tendon "
      "(compression > 0) [Pa], the initial tendon stress sigma_pi, the guaranteed strength fprg, Ep, the relative "
      "humidity rh [%] and the final pure relaxation [Pa]: elastic 3 sigma_b (n / 2 sigma_b), shrinkage eps_r Ep "
      "(4e-4 below 65 % of humidity, 3e-4 above; eps_r imposes it), creep 2 n sigma_b, viscoelastic total = "
      "shrinkage + creep + relaxation (1 - (shrinkage + creep) / (sigma_pi - 0.55 fprg)), the relaxation dropped "
      "once the stress loss exceeds that margin. Returns dict(elastic, shrinkage, creep, relaxation, viscoelastic, "
      "total = elastic + viscoelastic) [Pa].");
  mod.def(
      "bpel_lump_sum",
      [lump_dict, lump_options](double sigma_b, double sigma_pi, double fprg, double Ep, double rh, double rho1000,
                                double n, std::optional<double> eps_r, std::optional<double> mu0) {
        const LumpSumOptions o = lump_options(n, eps_r, mu0);
        return lump_dict(engine([&] { return bpel_lump_sum(sigma_b, sigma_pi, fprg, Ep, rh, rho1000, o); }));
      },
      py::arg("sigma_b"), py::arg("sigma_pi"), py::arg("fprg"), py::arg("Ep"), py::arg("rh"), py::arg("rho1000"),
      py::arg("n") = 6.0, py::arg("eps_r") = py::none(), py::arg("mu0") = py::none(),
      "Lump-sum deferred losses of the BPEL rules (arguments of ip2_lump_sum; rho1000 [%] in place of the final "
      "relaxation): elastic 3 sigma_b, shrinkage eps_r Ep (3e-4 below 65 % of humidity, 2e-4 above), creep 2.5 n "
      "sigma_b, relaxation 6/100 rho1000 (sigma_pi / fprg - mu0) sigma_pi (mu0 = bpel_mu0(rho1000) unless given, "
      "zero below it), viscoelastic total = shrinkage + creep + 5/6 relaxation. Returns the dict of ip2_lump_sum.");
  mod.def("bpel_mu0", [](double rho1000) { return engine([&] { return bpel_mu0(rho1000); }); }, py::arg("rho1000"),
          "mu0 of the BPEL relaxation by the band of rho1000 [%]: 0.43 up to 2.5, 0.35 up to 8, 0.30 above.");

  mod.def(
      "volterra_strain",
      [](const py::object& J, const std::vector<double>& t, const std::vector<double>& sigma) {
        const ComplianceFn f = compliance_of(J);
        return engine([&] { return volterra_strain(f, t, sigma); });
      },
      py::arg("J"), py::arg("t"), py::arg("sigma"),
      "Strain history of a stress history by the hereditary integral (trapezoidal); J = a concrete_law or a "
      "callable J(t, t0) [1/Pa]; sigma[0] applied at t[0] > 0.");
  mod.def(
      "volterra_stress",
      [](const py::object& J, const std::vector<double>& t, const std::vector<double>& eps) {
        const ComplianceFn f = compliance_of(J);
        return engine([&] { return volterra_stress(f, t, eps); });
      },
      py::arg("J"), py::arg("t"), py::arg("eps"), "Stress history of a strain history (step-by-step inversion).");
  mod.def(
      "relaxation_function",
      [](const py::object& J, double t0, const std::vector<double>& t) {
        const ComplianceFn f = compliance_of(J);
        return engine([&] { return relaxation_function(f, t0, t); });
      },
      py::arg("J"), py::arg("t0"), py::arg("t"), "R(t, t0): stress under a unit strain imposed at t0 (t[0] = t0).");
  mod.def("log_time_grid", &log_time_grid, py::arg("t0"), py::arg("dt_min"), py::arg("dt_max"), py::arg("n"));
  mod.def(
      "fit_kelvin_chain",
      [](const py::object& J, double t0, int n_units, double tau_min, double decades, double fit_min, double fit_max,
         int n_samples) {
        const ComplianceFn f = compliance_of(J);
        KelvinFitOptions o;
        o.n_units = n_units;
        o.tau_min = tau_min;
        o.decades = decades;
        o.fit_min = fit_min;
        o.fit_max = fit_max;
        o.n_samples = n_samples;
        const KelvinChain c = engine([&] { return fit_kelvin_chain(f, t0, o); });
        py::dict d;
        d["t0"] = c.t0;
        d["a0"] = c.a0;
        d["tau"] = c.tau;
        d["a"] = c.a;
        d["max_rel_error"] = c.fit_max_rel_error;
        return d;
      },
      py::arg("J"), py::arg("t0"), py::arg("n_units") = 8, py::arg("tau_min") = 1e-2, py::arg("decades") = 1.0,
      py::arg("fit_min") = 1e-2, py::arg("fit_max") = 1e5, py::arg("n_samples") = 120,
      "Non-negative Kelvin chain fitted to J(t0 + dt, t0): dict(a0 = J(t0, t0) exactly, tau, a, max_rel_error); "
      "J_chain(dt) = a0 + sum a_i (1 - exp(-dt / tau_i)).");
}
