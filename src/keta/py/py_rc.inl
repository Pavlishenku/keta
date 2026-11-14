// keta Python bindings — reinforced-concrete section analysis to EN 1992-1-1:
// SLS states I and II, cracking moment, crack width, tension
// stiffening and deflection, shrinkage curvature, ULS bending resistance and
// the N-M interaction curve; shear, torsion and punching resistances.
// Raw numbers only, never a verdict. Included by keta_py.cpp inside its
// anonymous namespace.
//
//   s = keta.RcSection(b=0.3, h=0.6, d=0.55, As=12.57e-4, fck=30e6)
//   s.cracked(N=0, M=150e3, phi=2.0)["sigma_s"], s.crack_width(0, 150e3)["wk"]
//   s.mrd(N=-500e3)["M_Rd"]

struct PyRc {
  RcSection s;
  RcMaterials m;
};

void bind_rc(py::module_& mod) {
  py::class_<PyRc>(mod, "RcSection",
                   "Rectangular or T reinforced-concrete section, EN 1992-1-1. z downward from the top fibre; N "
                   "tension-positive at the gross centroid; M > 0 sagging; stresses tension-positive.")
      .def(py::init([](double b, double h, double d, double As, double fck, double bf, double hf, double d2,
                       double As2, double fctm, double Ecm, double gamma_c, double alpha_cc, double fyk, double Es,
                       double gamma_s) {
             PyRc r;
             r.s.b = b;
             r.s.h = h;
             r.s.bf = bf;
             r.s.hf = hf;
             r.s.d = d;
             r.s.As = As;
             r.s.d2 = d2;
             r.s.As2 = As2;
             r.m.c.fck = fck;
             r.m.c.fctm = fctm;
             r.m.c.Ecm = Ecm;
             r.m.c.gamma_c = gamma_c;
             r.m.c.alpha_cc = alpha_cc;
             r.m.s.fyk = fyk;
             r.m.s.Es = Es;
             r.m.s.gamma_s = gamma_s;
             engine([&] { r.s.validate(); });
             return r;
           }),
           py::kw_only(), py::arg("b"), py::arg("h"), py::arg("d"), py::arg("As"), py::arg("fck"), py::arg("bf") = 0.0,
           py::arg("hf") = 0.0, py::arg("d2") = 0.0, py::arg("As2") = 0.0, py::arg("fctm") = 0.0, py::arg("Ecm") = 0.0,
           py::arg("gamma_c") = 1.5, py::arg("alpha_cc") = 1.0, py::arg("fyk") = 500e6, py::arg("Es") = 200e9,
           py::arg("gamma_s") = 1.15,
           "b web width, h depth, d / As tension steel, d2 / As2 compression steel, bf x hf top flange [m, m^2]; "
           "fck, fctm (0 = Table 3.1), Ecm (0 = EN value), fyk, Es [Pa]; partial factors.")
      .def(
          "uncracked",
          [](const PyRc& r, double N, double M, double phi) {
            const RcStateI q = engine([&] { return rc_uncracked(r.s, r.m, N, M, phi); });
            py::dict d;
            d["alpha_e"] = q.alpha_e;
            d["A"] = q.A;
            d["x"] = q.x;
            d["I"] = q.I;
            d["sigma_c_top"] = q.sigma_c_top;
            d["sigma_c_bot"] = q.sigma_c_bot;
            d["sigma_s"] = q.sigma_s;
            d["sigma_s2"] = q.sigma_s2;
            d["kappa"] = q.kappa;
            return d;
          },
          py::arg("N"), py::arg("M"), py::arg("phi") = 0.0,
          "SLS state I (transformed uncracked section, alpha_e = Es (1 + phi) / Ecm).")
      .def(
          "cracking_moment",
          [](const PyRc& r, double N, double phi, double fct_eff) {
            return engine([&] { return rc_cracking_moment(r.s, r.m, N, phi, fct_eff); });
          },
          py::arg("N") = 0.0, py::arg("phi") = 0.0, py::arg("fct_eff") = 0.0,
          "Sagging moment bringing the bottom fibre of the state-I section to fct_eff (0 = fctm).")
      .def(
          "cracked",
          [](const PyRc& r, double N, double M, double phi) {
            const RcStateII q = engine([&] { return rc_cracked(r.s, r.m, N, M, phi); });
            py::dict d;
            d["alpha_e"] = q.alpha_e;
            d["x"] = q.x;
            d["I"] = q.I;
            d["sigma_c_top"] = q.sigma_c_top;
            d["sigma_s"] = q.sigma_s;
            d["sigma_s2"] = q.sigma_s2;
            d["kappa"] = q.kappa;
            d["lever_arm"] = q.lever_arm;
            return d;
          },
          py::arg("N"), py::arg("M"), py::arg("phi") = 0.0,
          "SLS state II (concrete in tension ignored, linear concrete in compression).")
      .def(
          "crack_width",
          [](const PyRc& r, double N, double M, double phi, double cover, double bar_dia, double spacing, double k1,
             double k2, double k3, double k4, double kt, double fct_eff) {
            RcCrackParams p;
            p.cover = cover;
            p.bar_dia = bar_dia;
            p.spacing = spacing;
            p.k1 = k1;
            p.k2 = k2;
            p.k3 = k3;
            p.k4 = k4;
            p.kt = kt;
            p.fct_eff = fct_eff;
            const RcCrackResult q = engine([&] { return rc_crack_width(r.s, r.m, p, N, M, phi); });
            py::dict d;
            d["x"] = q.x;
            d["sigma_s"] = q.sigma_s;
            d["hc_ef"] = q.hc_ef;
            d["Ac_eff"] = q.Ac_eff;
            d["rho_p_eff"] = q.rho_p_eff;
            d["d_eps"] = q.d_eps;
            d["s_r_max"] = q.s_r_max;
            d["wide_spacing"] = q.wide_spacing;
            d["wk"] = q.wk;
            return d;
          },
          py::arg("N"), py::arg("M"), py::arg("phi") = 0.0, py::arg("cover") = 0.030, py::arg("bar_dia") = 0.016,
          py::arg("spacing") = 0.0, py::arg("k1") = 0.8, py::arg("k2") = 0.5, py::arg("k3") = 3.4,
          py::arg("k4") = 0.425, py::arg("kt") = 0.6, py::arg("fct_eff") = 0.0,
          "Crack width §7.3.4: wk = s_r,max (eps_sm - eps_cm), eq. 7.8, 7.9, 7.11 (7.14 when the spacing exceeds "
          "5 (c + phi/2)).")
      .def(
          "curvature",
          [](const PyRc& r, double N, double M, double phi, double beta) {
            const RcCurvature q = engine([&] { return rc_curvature(r.s, r.m, N, M, phi, beta); });
            py::dict d;
            d["M_cr"] = q.M_cr;
            d["kappa_I"] = q.kappa_I;
            d["kappa_II"] = q.kappa_II;
            d["zeta"] = q.zeta;
            d["kappa"] = q.kappa;
            return d;
          },
          py::arg("N"), py::arg("M"), py::arg("phi") = 0.0, py::arg("beta") = 0.5,
          "Tension-stiffening curvature §7.4.3 eq. 7.18 / 7.19 (beta 1 short-term single load, 0.5 sustained).")
      .def(
          "deflection_ss",
          [](const PyRc& r, double L, const std::function<double(double)>& M_of_x, double N, double phi, double beta,
             int n, double eps_cs, double phi_cr) {
            return engine([&] {
              return eps_cs != 0.0 ? rc_deflection_ss_cs(r.s, r.m, L, M_of_x, phi, eps_cs, beta, phi_cr)
                                   : rc_deflection_ss(r.s, r.m, L, M_of_x, N, phi, beta, n);
            });
          },
          py::arg("L"), py::arg("M_of_x"), py::arg("N") = 0.0, py::arg("phi") = 0.0, py::arg("beta") = 0.5,
          py::arg("n_per_half") = 40, py::arg("eps_cs") = 0.0, py::arg("phi_cr") = -1.0,
          "Mid-span deflection of a simply supported beam under M(x) with the eq. 7.18 curvature (unit-load "
          "integral); eps_cs != 0 adds the shrinkage curvature, interpolated by the same zeta (N = 0 then).")
      .def(
          "shrinkage_curvature",
          [](const PyRc& r, double eps_cs, double phi, bool cracked) {
            return engine([&] { return rc_shrinkage_curvature(r.s, r.m, eps_cs, phi, cracked); });
          },
          py::arg("eps_cs"), py::arg("phi") = 0.0, py::arg("cracked") = false,
          "Eq. 7.21: 1/r_cs = eps_cs alpha_e S / I (state I or II), sagging-positive.")
      .def(
          "mrd",
          [](const PyRc& r, double N) {
            const RcUls q = engine([&] { return rc_mrd(r.s, r.m, N); });
            py::dict d;
            d["fcd"] = q.fcd;
            d["fyd"] = q.fyd;
            d["lambda"] = q.lambda;
            d["eta"] = q.eta;
            d["eps_cu3"] = q.eps_cu3;
            d["x"] = q.x;
            d["block"] = q.block;
            d["eps_s"] = q.eps_s;
            d["sigma_s"] = q.sigma_s;
            d["sigma_s2"] = q.sigma_s2;
            d["Fc"] = q.Fc;
            d["M_Rd"] = q.M_Rd;
            d["steel_yields"] = q.steel_yields;
            return d;
          },
          py::arg("N") = 0.0,
          "ULS sagging resistance with the rectangular block §3.1.7(3) at the axial force N (tension +), about the "
          "gross centroid.")
      .def(
          "interaction",
          [](const PyRc& r, int n) {
            if (n < 2) throw py::value_error("interaction: n >= 2 points");
            // bounds: full tension of the steel, full compression of the block
            const RcUls u = engine([&] { return rc_mrd(r.s, r.m, 0.0); });
            const double Nt = (r.s.As + r.s.As2) * u.fyd;
            const double Nc = -(u.eta * u.fcd * r.s.area() + (r.s.As + r.s.As2) * u.fyd);
            std::vector<std::pair<double, double>> out;
            for (int i = 0; i < n; ++i) {
              const double N = Nc + (Nt - Nc) * (i + 0.5) / n;
              try {
                out.emplace_back(N, rc_mrd(r.s, r.m, N).M_Rd);
              } catch (const std::exception&) {
                // beyond the capacity the section has no equilibrium: skipped
              }
            }
            return out;
          },
          py::arg("n") = 41,
          "Sagging branch of the N-M interaction curve: (N, M_Rd) at n axial forces spread between the squash "
          "and the tension capacities (points with no equilibrium are skipped).");

  // ---- shear, torsion, punching (raw resistances, ec2_design.h) ----
  mod.def(
      "vrd_c",
      [](double bw, double d, double Asl, double fck, double gamma_c, double sigma_cp, double fcd, bool din) {
        ec2::VrdcInput in;
        in.bw = bw;
        in.d = d;
        in.Asl = Asl;
        in.fck = fck;
        in.gamma_c = gamma_c;
        in.sigma_cp = sigma_cp;
        in.fcd = fcd;
        return engine([&] { return din ? ec2::vrd_c_din(in) : ec2::vrd_c_en(in); });
      },
      py::kw_only(), py::arg("bw"), py::arg("d"), py::arg("Asl"), py::arg("fck"), py::arg("gamma_c") = 1.5,
      py::arg("sigma_cp") = 0.0, py::arg("fcd") = 0.0, py::arg("din") = false,
      "§6.2.2 eq. 6.2a/6.2b shear resistance without shear reinforcement [N] (din: the German NA).");
  mod.def("asw_s", &ec2::asw_s, py::arg("VEd"), py::arg("z"), py::arg("fywd"), py::arg("cot_theta"),
          "Eq. 6.8 required Asw / s = VEd / (z fywd cot(theta)) [m^2/m].");
  mod.def("vrd_max", &ec2::vrd_max, py::arg("bw"), py::arg("z"), py::arg("nu1"), py::arg("fcd"), py::arg("cot_theta"),
          py::arg("alpha_cw") = 1.0, "Eq. 6.9 strut crushing VRd,max = alpha_cw bw z nu1 fcd / (cot + tan) [N].");
  mod.def("nu_shear", &ec2::nu_en, py::arg("fck"), "Eq. 6.6N nu = 0.6 (1 - fck / 250 MPa).");
  mod.def("asw_min_s", &ec2::asw_min_s_en, py::arg("fck"), py::arg("fyk"), py::arg("bw"), py::arg("alpha") = kPi / 2,
          "Eq. 9.4/9.5N minimum Asw / s = 0.08 sqrt(fck) / fyk bw sin(alpha) [m^2/m].");
  mod.def(
      "torsion_reinforcement",
      [](double TEd, double b, double h, double tef, double fyd, double cot_theta, double nu, double fcd) {
        const ec2::TorsionBox t = ec2::torsion_box(b, h, tef);
        py::dict d;
        d["Ak"] = t.Ak;
        d["uk"] = t.uk;
        d["asw_s"] = ec2::asw_s_torsion(TEd, t, fyd, cot_theta);
        d["asl"] = ec2::asl_torsion(TEd, t, fyd, cot_theta);
        d["TRd_max"] = nu > 0 ? ec2::trd_max(nu, fcd, t, cot_theta) : 0.0;
        return d;
      },
      py::kw_only(), py::arg("TEd"), py::arg("b"), py::arg("h"), py::arg("tef"), py::arg("fyd"),
      py::arg("cot_theta") = 1.0, py::arg("nu") = 0.0, py::arg("fcd") = 0.0,
      "§6.3.2 thin-walled box: Ak, uk, stirrups Asw/s, longitudinal Asl (eq. 6.28) and TRd,max (eq. 6.30).");
  mod.def(
      "punching",
      [](double VEd, double a, double b, double d, double rho_ly, double rho_lz, double fck, double beta,
         double gamma_c, double sigma_cp) {
        ec2::PunchingInput in;
        in.d = d;
        in.rho_ly = rho_ly;
        in.rho_lz = rho_lz;
        in.fck = fck;
        in.gamma_c = gamma_c;
        in.sigma_cp = sigma_cp;
        in.u0 = 2 * (a + b);
        py::dict o;
        const double u1 = ec2::punching_u1_rect(a, b, d);
        o["u1"] = u1;
        o["vRd_c"] = engine([&] { return ec2::vrd_c_punching_en(in); });
        o["vEd"] = ec2::ved_punching(beta, VEd, u1, d);
        return o;
      },
      py::kw_only(), py::arg("VEd"), py::arg("a"), py::arg("b"), py::arg("d"), py::arg("rho_ly"), py::arg("rho_lz"),
      py::arg("fck"), py::arg("beta") = 1.0, py::arg("gamma_c") = 1.5, py::arg("sigma_cp") = 0.0,
      "§6.4 punching of an interior a x b column: u1 = 2 (a + b) + 4 pi d, vRd,c (eq. 6.47), vEd = beta VEd / "
      "(u1 d) [Pa].");

  // ---- prestressed sections (ec2_section.h and fiber.h) ----
  // a memoryless law from {'kind': ..., params}: 'elastic' (E), 'concrete_pr'
  // (fcd, eps_c2, eps_cu2, n), 'steel' (E, fy, eps_u) or a design law of
  // design_law_of ('strand', 'sargin', 'steel_hardening', ...)
  auto law_of = [](const py::dict& d, const char* what) {
    if (!d.contains("kind")) throw py::value_error(std::string(what) + ": the law dict needs a 'kind'");
    const std::string k = d["kind"].cast<std::string>();
    auto g = [&](const char* key) {
      if (!d.contains(key)) throw py::value_error(std::string(what) + ": law '" + k + "' needs '" + key + "'");
      return d[key].cast<double>();
    };
    auto has = [&](const char* key) { return d.contains(key); };
    auto opt = [&](const char* key, double v) { return has(key) ? g(key) : v; };
    return engine([&]() -> UniaxialLaw {
      if (k == "elastic") return UniaxialLaw::elastic(g("E"));
      if (k == "concrete_pr") return UniaxialLaw::concrete_pr(g("fcd"), opt("eps_c2", 2e-3), opt("eps_cu2", 3.5e-3), opt("n", 2.0));
      if (k == "steel") return UniaxialLaw::steel(g("E"), g("fy"), opt("eps_u", 0.0));
      if (auto l = design_law_of(k, g, has)) return *l;
      throw py::value_error(std::string(what) + ": law kind 'elastic', 'concrete_pr', 'steel', 'strand', "
                            "'steel_hardening', 'sargin', 'concrete_tension_softening' or 'bfup_t3'");
    });
  };
  // the exact section with its bars and tendons
  auto poly_section = [law_of](const std::vector<std::vector<std::pair<double, double>>>& rings,
                               const std::vector<std::vector<double>>& bars, double fck, double fcd,
                               const py::list& tendons, const py::object& strand, double fyd, double ftd,
                               double eps_ud, double Es) {
    ec2::RcPolySection s;
    for (const auto& r : rings) {
      std::vector<Vec2> ring;
      for (const auto& [y, z] : r) ring.emplace_back(y, z);
      s.rings.push_back(ring);
    }
    // a bar is (y, z, weight) or (y, z, weight, eps0), eps0 its initial strain
    for (const std::vector<double>& b : bars) {
      if (b.size() != 3 && b.size() != 4) throw py::value_error("bars: (y, z, weight) or (y, z, weight, eps0)");
      s.bars.push_back({b[0], b[1], b[2], b.size() == 4 ? b[3] : 0.0});
    }
    s.fck = fck;
    s.fcd = fcd;
    s.steel.fyd = fyd;
    s.steel.ftd = ftd;
    s.steel.eps_ud = eps_ud;
    s.steel.Es = Es;
    if (!strand.is_none()) s.strand = law_of(strand.cast<py::dict>(), "prestressed section");
    for (const py::handle& h : tendons) {
      const py::dict t = py::reinterpret_borrow<py::dict>(h);
      ec2::RcPolySection::Tendon tt;
      tt.y = t["y"].cast<double>();
      tt.z = t["z"].cast<double>();
      tt.Ap = t["Ap"].cast<double>();
      tt.eps_p0 = t["eps_p0"].cast<double>();
      if (t.contains("bonded")) tt.bonded = t["bonded"].cast<bool>();
      if (t.contains("dsigma_uls")) tt.dsigma_uls = t["dsigma_uls"].cast<double>();
      if (t.contains("displaced_area")) tt.displaced_area = t["displaced_area"].cast<double>();
      if (t.contains("displaced_dy")) tt.displaced_dy = t["displaced_dy"].cast<double>();
      if (t.contains("displaced_dz")) tt.displaced_dz = t["displaced_dz"].cast<double>();
      s.tendons.push_back(tt);
    }
    return s;
  };
  using Rings = std::vector<std::vector<std::pair<double, double>>>;
  using Bars = std::vector<std::vector<double>>;
  mod.def(
      "prestressed_section_resistance",
      [poly_section, law_of](const Rings& rings, double fck, double fcd, const py::list& tendons, const Bars& bars,
                             double omega, double N, double theta, const py::object& strand, double fyd,
                             double ftd, double eps_ud, double Es, const py::object& concrete,
                             bool deduct_displaced) {
        ec2::RcPolySection s = poly_section(rings, bars, fck, fcd, tendons, strand, fyd, ftd, eps_ud, Es);
        if (!concrete.is_none()) s.concrete = law_of(concrete.cast<py::dict>(), "prestressed section");
        s.deduct_displaced = deduct_displaced;
        const ec2::SectionDesign r = engine([&] { return ec2::moment_resistance(s, omega, N, theta); });
        py::dict d;
        d["feasible"] = r.feasible;
        d["M_R"] = std::hypot(r.R(1), r.R(2));
        d["resultants"] = std::make_tuple(r.R(0), r.R(1), r.R(2));
        d["plane"] = std::make_tuple(r.plane.e0, r.plane.ky, r.plane.kz);
        d["tau"] = r.tau;
        std::vector<double> sp;
        if (r.feasible)
          for (const auto& t : s.tendons) sp.push_back(engine([&] { return ec2::tendon_stress(s, t, r.plane); }));
        d["tendon_stress"] = sp;
        return d;
      },
      py::arg("rings"), py::arg("fck"), py::arg("fcd"), py::arg("tendons"), py::arg("bars") = Bars{},
      py::arg("omega") = 1.0, py::arg("N") = 0.0, py::arg("theta") = 0.0, py::arg("strand") = py::none(),
      py::arg("fyd") = 500e6 / 1.15, py::arg("ftd") = 525e6 / 1.15, py::arg("eps_ud") = 0.025, py::arg("Es") = 200e9,
      py::arg("concrete") = py::none(), py::arg("deduct_displaced") = false,
      "EN 1992-1-1 resisting moment of a prestressed polygonal section (rings in (y, z): outer CCW, holes CW) at "
      "N, along the direction theta of (My, Mz). tendons: [dict(y, z, Ap, eps_p0, bonded=True, dsigma_uls=100e6)]: a "
      "bonded tendon follows the plane from its initial strain eps_p0 (Fig. 6.1), an unbonded one keeps its "
      "stress plus dsigma_uls (§5.10.8(2)); bars [(y, z, weight)] of area omega weight; strand: the tendon law "
      "(dict kind='strand', ...; default Y1860S7). concrete: law dict of the concrete (kind='concrete_pr' with "
      "fcd, eps_c2, eps_cu2, n — n = 1 is the bilinear relation, see bilinear_concrete — or kind='sargin'); None: "
      "the parabola-rectangle of fck at fcd. deduct_displaced: the concrete displaced by the tendons "
      "(displaced_area) and by the bars (their own area) is deducted. Returns feasible, M_R, resultants (N, My, "
      "Mz), plane (e0, ky, kz), tau (pivot A [0, 1], B, C) and tendon_stress (keta/ec2_section.h).");
  mod.def(
      "prestressed_section_diagram",
      [poly_section, law_of](const Rings& rings, double fck, double fcd, const py::list& tendons, double alpha,
                             int n_points, const Bars& bars, double omega, const py::object& strand, double fyd,
                             double ftd, double eps_ud, double Es, const py::object& concrete,
                             bool deduct_displaced) {
        ec2::RcPolySection s = poly_section(rings, bars, fck, fcd, tendons, strand, fyd, ftd, eps_ud, Es);
        if (!concrete.is_none()) s.concrete = law_of(concrete.cast<py::dict>(), "prestressed section");
        s.deduct_displaced = deduct_displaced;
        std::vector<Vec3T> out;
        for (const Vec3& v : engine([&] { return ec2::interaction_diagram(s, omega, alpha, n_points); }))
          out.emplace_back(v(0), v(1), v(2));
        return out;
      },
      py::arg("rings"), py::arg("fck"), py::arg("fcd"), py::arg("tendons"), py::arg("alpha"), py::arg("n_points") = 61,
      py::arg("bars") = Bars{}, py::arg("omega") = 1.0, py::arg("strand") = py::none(), py::arg("fyd") = 500e6 / 1.15,
      py::arg("ftd") = 525e6 / 1.15, py::arg("eps_ud") = 0.025, py::arg("Es") = 200e9,
      py::arg("concrete") = py::none(), py::arg("deduct_displaced") = false,
      "N-M interaction diagram of a prestressed section, compression on the side of (cos alpha, sin alpha): "
      "[(N, My, Mz)] of the ultimate planes from pivot A to the squash, N decreasing (same section arguments as "
      "prestressed_section_resistance, concrete law and displaced concrete included).");
  mod.def(
      "bilinear_concrete",
      [](double fck, double fcd) {
        const UniaxialLaw l = engine([&] { return ec2::bilinear_concrete(fck, fcd); });
        py::dict d;
        d["kind"] = "concrete_pr";
        d["fcd"] = l.fcd;
        d["eps_c2"] = l.eps_c2;
        d["eps_cu2"] = l.eps_cu2;
        d["n"] = l.n_exp;
        return d;
      },
      py::arg("fck"), py::arg("fcd"),
      "The bilinear stress-strain relation of EN 1992-1-1 §3.1.7(2) as a law dict: a line to fcd at eps_c3, a "
      "plateau to eps_cu3 (Table 3.1 for fck) — the parabola-rectangle of exponent 1, dict(kind='concrete_pr', "
      "fcd, eps_c2=eps_c3, eps_cu2=eps_cu3, n=1).");
  mod.def(
      "prestress_decomposition",
      [poly_section](const py::list& tendons, const Vec3T& total, const py::object& strand) {
        const ec2::RcPolySection s = poly_section({}, {}, 30e6, 20e6, tendons, strand, 500e6 / 1.15, 525e6 / 1.15, 0.025, 200e9);
        const ec2::PrestressDecomposition d = engine([&] {
          return ec2::decompose_prestress(s, Vec3(std::get<0>(total), std::get<1>(total), std::get<2>(total)));
        });
        py::dict o;
        o["isostatic"] = Vec3T(d.isostatic(0), d.isostatic(1), d.isostatic(2));
        o["hyperstatic"] = Vec3T(d.hyperstatic(0), d.hyperstatic(1), d.hyperstatic(2));
        return o;
      },
      py::arg("tendons"), py::arg("total") = Vec3T(0.0, 0.0, 0.0), py::arg("strand") = py::none(),
      "Isostatic action of the tendons on the section, (N, My, Mz) = sum (-P, -P z, +P y) with P = Ap "
      "sigma(eps_p0), and the hyperstatic remainder total - isostatic of a total prestress effect (N, My, Mz) "
      "from a frame analysis.");
  mod.def(
      "tendon_initial_strain",
      [law_of](double sigma_pm, double sigma_cp, double Ecm, const py::object& strand) {
        const UniaxialLaw l = strand.is_none() ? ec2::RcPolySection{}.strand : law_of(strand.cast<py::dict>(), "tendon_initial_strain");
        return engine([&] { return ec2::tendon_initial_strain(l, sigma_pm, sigma_cp, Ecm); });
      },
      py::arg("sigma_pm"), py::arg("sigma_cp") = 0.0, py::arg("Ecm") = 0.0, py::arg("strand") = py::none(),
      "Initial strain Delta eps_p(0) of a tendon at stress sigma_pm next to concrete at sigma_cp (< 0 in "
      "compression): the strain of sigma_pm on the law plus -sigma_cp / Ecm.");
  // The prestressed fiber rectangle of the Python calls: rect = (b, h, ny, nz), a concrete law dict,
  // tendons [dict(y, z, A, eps_p0, bonded=True, law=strand dict)] and rebars [dict(y, z, A, law, eps0=0)].
  const auto fiber_rectangle = [law_of](const std::tuple<double, double, int, int>& rect, const py::dict& concrete,
                                        const py::list& tendons, const py::list& rebars, const char* who) {
    FiberSection f = engine([&] {
      return FiberSection::rect(std::get<0>(rect), std::get<1>(rect), std::get<2>(rect), std::get<3>(rect),
                                law_of(concrete, who));
    });
    for (const py::handle& h : rebars) {
      const py::dict r = py::reinterpret_borrow<py::dict>(h);
      const UniaxialLaw l = law_of(r["law"].cast<py::dict>(), who);
      engine([&] {
        f.add_rebar(r["y"].cast<double>(), r["z"].cast<double>(), r["A"].cast<double>(), l,
                    r.contains("eps0") ? r["eps0"].cast<double>() : 0.0);
        return 0;
      });
    }
    for (const py::handle& h : tendons) {
      const py::dict t = py::reinterpret_borrow<py::dict>(h);
      const UniaxialLaw l = t.contains("law") ? law_of(t["law"].cast<py::dict>(), who) : ec2::RcPolySection{}.strand;
      engine([&] {
        f.add_tendon(t["y"].cast<double>(), t["z"].cast<double>(), t["A"].cast<double>(), l,
                     t["eps_p0"].cast<double>(), t.contains("bonded") ? t["bonded"].cast<bool>() : true);
        return 0;
      });
    }
    return f;
  };
  mod.def(
      "prestressed_m_kappa",
      [fiber_rectangle](const std::tuple<double, double, int, int>& rect, const py::dict& concrete,
                        const py::list& tendons, const std::vector<double>& kappas, double N, const py::list& rebars) {
        const FiberSection f = fiber_rectangle(rect, concrete, tendons, rebars, "prestressed_m_kappa");
        std::vector<Vec3T> out;
        for (const auto& p : engine([&] { return f.m_kappa(N, kappas); })) out.emplace_back(p.kappa, p.moment, p.eps0);
        return out;
      },
      py::arg("rect"), py::arg("concrete"), py::arg("tendons"), py::arg("kappas"), py::arg("N") = 0.0,
      py::arg("rebars") = py::list(),
      "Moment-curvature of a prestressed fiber rectangle rect = (b, h, ny, nz) at constant N: concrete law dict, "
      "tendons [dict(y, z, A, eps_p0, bonded=True, law=strand dict)], rebars [dict(y, z, A, law, eps0=0)]; returns "
      "[(kappa_y, My, eps0)] (keta/fiber.h; unbonded tendons keep their initial force).");
  mod.def(
      "prestressed_ultimate_point",
      [fiber_rectangle](const std::tuple<double, double, int, int>& rect, const py::dict& concrete,
                        const py::list& tendons, double N, const py::list& rebars, double alpha,
                        const py::object& eps_steel_max) {
        const FiberSection f = fiber_rectangle(rect, concrete, tendons, rebars, "prestressed_ultimate_point");
        const double steel_max =
            eps_steel_max.is_none() ? std::numeric_limits<double>::infinity() : eps_steel_max.cast<double>();
        const UltimatePoint p = engine([&] { return ultimate_point(f, N, alpha, steel_max); });
        py::dict o;
        o["kappa_y"] = p.kappa * std::cos(alpha);
        o["kappa_z"] = p.kappa * std::sin(alpha);
        o["My"] = p.my;
        o["Mz"] = p.mz;
        o["eps0"] = p.eps0;
        return o;
      },
      py::arg("rect"), py::arg("concrete"), py::arg("tendons"), py::arg("N") = 0.0, py::arg("rebars") = py::list(),
      py::arg("alpha") = 0.0, py::arg("eps_steel_max") = py::none(),
      "Ultimate point of the moment-curvature of a prestressed fiber rectangle (arguments of prestressed_m_kappa) "
      "at constant N, bending direction alpha (kappa_y = k cos alpha, kappa_z = k sin alpha): the largest curvature "
      "that keeps every fibre and bonded tendon within the strain bounds of its law. eps_steel_max also bounds the "
      "bilinear steel fibres (None: their own bound only). Returns dict(kappa_y, kappa_z, My, Mz, eps0).");

  mod.def(
      "fiber_resultants",
      [fiber_rectangle](const std::tuple<double, double, int, int>& rect, const py::dict& concrete, const Vec3T& plane,
                        const py::list& tendons, const py::list& rebars) {
        const FiberSection f = fiber_rectangle(rect, concrete, tendons, rebars, "fiber_resultants");
        const Vec3 r = engine([&] { return f.resultants(Vec3(std::get<0>(plane), std::get<1>(plane), std::get<2>(plane))); });
        return Vec3T(r(0), r(1), r(2));
      },
      py::arg("rect"), py::arg("concrete"), py::arg("plane"), py::arg("tendons") = py::list(),
      py::arg("rebars") = py::list(),
      "Resultants (N, My, Mz) of a fiber rectangle (arguments of prestressed_m_kappa) under the strain plane "
      "(eps0, kappa_y, kappa_z): eps(y, z) = eps0 + kappa_y z - kappa_z y, each fibre, bar and bonded tendon on "
      "its law.");
  mod.def(
      "fiber_state",
      [fiber_rectangle](const std::tuple<double, double, int, int>& rect, const py::dict& concrete, const Vec3T& action,
                        const py::list& tendons, const py::list& rebars) {
        const FiberSection f = fiber_rectangle(rect, concrete, tendons, rebars, "fiber_state");
        const Vec3 e = engine([&] { return f.solve_state(Vec3(std::get<0>(action), std::get<1>(action), std::get<2>(action))); });
        return Vec3T(e(0), e(1), e(2));
      },
      py::arg("rect"), py::arg("concrete"), py::arg("action"), py::arg("tendons") = py::list(),
      py::arg("rebars") = py::list(),
      "Strain plane (eps0, kappa_y, kappa_z) of a fiber rectangle (arguments of prestressed_m_kappa) in "
      "equilibrium with the action (N, My, Mz), every fibre within the strain bounds of its law; an action out of "
      "the capacity of the section is refused.");
  mod.def(
      "fiber_interaction",
      [fiber_rectangle](const std::tuple<double, double, int, int>& rect, const py::dict& concrete, double N,
                        int n_alpha, double eps_steel_max, const py::list& tendons, const py::list& rebars) {
        if (n_alpha < 1) throw py::value_error("fiber_interaction: n_alpha >= 1");
        const FiberSection f = fiber_rectangle(rect, concrete, tendons, rebars, "fiber_interaction");
        std::vector<Vec3T> out;
        for (const InteractionPoint& p : engine([&] { return interaction_curve(f, N, n_alpha, eps_steel_max); }))
          out.emplace_back(p.alpha, p.my, p.mz);
        return out;
      },
      py::arg("rect"), py::arg("concrete"), py::arg("N") = 0.0, py::arg("n_alpha") = 16,
      py::arg("eps_steel_max") = 0.01, py::arg("tendons") = py::list(), py::arg("rebars") = py::list(),
      "My-Mz interaction curve of a fiber rectangle (arguments of prestressed_m_kappa) at constant N: the ultimate "
      "point of n_alpha bending directions spread over the turn, alpha = 2 pi i / n_alpha (kappa_y = k cos alpha, "
      "kappa_z = k sin alpha). eps_steel_max bounds the strain of the plateau steels. Returns [(alpha, My, Mz)].");

  // ---- prestress losses at the level of a section (keta/section_losses.h)
  // the concrete rings and the tendons, without reinforcement or design strengths
  const auto loss_section = [poly_section](const Rings& rings, const py::list& tendons, const py::object& strand) {
    return poly_section(rings, {}, 30e6, 20e6, tendons, strand, 500e6 / 1.15, 525e6 / 1.15, 0.025, 200e9);
  };
  const auto loss_dicts = [](const std::vector<ec2::TendonLoss>& losses) {
    py::list out;
    for (const ec2::TendonLoss& l : losses) {
      py::dict d;
      d["stress"] = l.stress;
      d["elastic_shortening"] = l.elastic_shortening;
      d["creep"] = l.creep;
      d["shrinkage"] = l.shrinkage;
      d["relaxation"] = l.relaxation;
      d["total"] = l.total;
      d["eps_p0_after"] = l.eps_p0_after;
      out.append(d);
    }
    return out;
  };
  // ---- service states of a section (keta/section_losses.h, keta/section_service.h)
  // the section of a service state: rings, tendons, bars [(y, z, weight)] of area omega weight
  const auto service_section = [poly_section](const Rings& rings, const py::list& tendons, const py::object& strand,
                                              const Bars& bars, double fy, double ft, double eps_ud, double Es) {
    return poly_section(rings, bars, 30e6, 20e6, tendons, strand, fy, ft, eps_ud, Es);
  };
  const auto passive = [](const Bars& bars, double omega, bool displace) {
    return ec2::PassiveSteel{bars.empty() ? 0.0 : omega, displace};
  };
  const auto points_of = [](const std::vector<std::pair<double, double>>& points) {
    std::vector<Vec2> at;
    for (const auto& [y, z] : points) at.emplace_back(y, z);
    return at;
  };
  const auto state_dict = [](const ec2::PrestressServiceState& st) {
    const auto fibre = [](const ec2::FibreStress& f) {
      py::dict d;
      d["y"] = f.y;
      d["z"] = f.z;
      d["stress"] = f.stress;
      return d;
    };
    py::dict d;
    d["tendon_stress"] = st.tendon_stress;
    d["concrete_stress"] = st.concrete_stress;
    d["point_stress"] = st.point_stress;
    d["fibre_max"] = fibre(st.fibre_max);
    d["fibre_min"] = fibre(st.fibre_min);
    d["bar_stress"] = st.bar_stress;
    d["plane"] = std::make_tuple(st.plane.e0, st.plane.ky, st.plane.kz);
    d["decompressed"] = st.decompressed;
    return d;
  };
  using Points = std::vector<std::pair<double, double>>;
  mod.def(
      "prestress_service_state",
      [service_section, passive, points_of, state_dict](const Rings& rings, const py::list& tendons, double Ecm,
                                                         const Vec3T& sustained, const py::object& strand,
                                                         const Points& points, const Bars& bars, double omega,
                                                         double Es, bool bars_displace_concrete) {
        const ec2::RcPolySection s = service_section(rings, tendons, strand, bars, 500e6, 525e6, 0.025, Es);
        return state_dict(engine([&] {
          return ec2::prestress_service_state(
              s, Ecm, Vec3(std::get<0>(sustained), std::get<1>(sustained), std::get<2>(sustained)), points_of(points),
              passive(bars, omega, bars_displace_concrete));
        }));
      },
      py::arg("rings"), py::arg("tendons"), py::arg("Ecm"), py::arg("sustained") = Vec3T(0.0, 0.0, 0.0),
      py::arg("strand") = py::none(), py::arg("points") = Points{}, py::arg("bars") = Bars{}, py::arg("omega") = 1.0,
      py::arg("Es") = 200e9, py::arg("bars_displace_concrete") = false,
      "Elastic service state of a prestressed or reinforced section, concrete uncracked (the concrete of the "
      "rings, modulus Ecm; tendons [dict(y, z, Ap, eps_p0, bonded=True, displaced_area=0, displaced_dy=0, "
      "displaced_dz=0)] of modulus Ep; bars [(y, z, weight)] of area omega weight and modulus Es) under its "
      "prestress and sustained = (N, My, Mz): dict(tendon_stress, concrete_stress, point_stress, fibre_max, "
      "fibre_min, bar_stress, plane, decompressed), tension positive. concrete_stress is the concrete stress at "
      "each tendon, point_stress at each of points = [(y, z)] (a point outside the concrete is refused); "
      "fibre_max and fibre_min = dict(y, z, stress) are the most tensioned and the most compressed vertex of the "
      "rings; plane = (e0, ky, kz) of the strain e0 + ky z - kz y. A bonded tendon carries Ep (eps_p0 + concrete "
      "strain), an unbonded one Ep eps_p0, a bar Es times the concrete strain. A tendon takes away the disc of "
      "concrete of area displaced_area centred (displaced_dy, displaced_dz) from it (an empty duct, or its own "
      "area once grouted); bars_displace_concrete does the same for each bar.");
  mod.def(
      "cracked_service_state",
      [service_section, passive, points_of, state_dict](const Rings& rings, const py::list& tendons, double Ecm,
                                                         const Vec3T& sustained, const py::object& strand,
                                                         const Points& points, const Bars& bars, double omega,
                                                         double Es, bool bars_displace_concrete) {
        const ec2::RcPolySection s = service_section(rings, tendons, strand, bars, 500e6, 525e6, 0.025, Es);
        return state_dict(engine([&] {
          return ec2::cracked_service_state(
              s, Ecm, Vec3(std::get<0>(sustained), std::get<1>(sustained), std::get<2>(sustained)), points_of(points),
              passive(bars, omega, bars_displace_concrete));
        }));
      },
      py::arg("rings"), py::arg("tendons"), py::arg("Ecm"), py::arg("sustained") = Vec3T(0.0, 0.0, 0.0),
      py::arg("strand") = py::none(), py::arg("points") = Points{}, py::arg("bars") = Bars{}, py::arg("omega") = 1.0,
      py::arg("Es") = 200e9, py::arg("bars_displace_concrete") = false,
      "Cracked service state of a prestressed or reinforced section (same arguments and result as "
      "prestress_service_state): the concrete carries Ecm times its strain where it is compressed and nothing "
      "where it is stretched, tendons and bars stay elastic. The strain plane is the one whose resultants are "
      "sustained = (N, My, Mz); fibre_max and fibre_min are the vertices of largest and smallest strain. A section "
      "left compressed everywhere is in its uncracked state. An action the section cannot carry is refused.");
  mod.def(
      "nonlinear_service_state",
      [service_section, passive, points_of, state_dict, law_of](
          const Rings& rings, const py::list& tendons, const py::dict& concrete, const Vec3T& sustained,
          const py::object& strand, const Points& points, const Bars& bars, double omega, double Es, double fy,
          double ft, double eps_ud, bool bars_displace_concrete) {
        ec2::RcPolySection s = service_section(rings, tendons, strand, bars, fy, ft, eps_ud, Es);
        // a service state reads the strand at its characteristic strengths
        if (strand.is_none()) s.strand = UniaxialLaw::strand(195e9, 1640e6, 1860e6, 0.035, 0.0315, true, 1.0);
        const UniaxialLaw law = law_of(concrete, "nonlinear service state");
        return state_dict(engine([&] {
          return ec2::nonlinear_service_state(
              s, law, Vec3(std::get<0>(sustained), std::get<1>(sustained), std::get<2>(sustained)), points_of(points),
              passive(bars, omega, bars_displace_concrete));
        }));
      },
      py::arg("rings"), py::arg("tendons"), py::arg("concrete"), py::arg("sustained") = Vec3T(0.0, 0.0, 0.0),
      py::arg("strand") = py::none(), py::arg("points") = Points{}, py::arg("bars") = Bars{}, py::arg("omega") = 1.0,
      py::arg("Es") = 200e9, py::arg("fy") = 500e6, py::arg("ft") = 525e6, py::arg("eps_ud") = 0.025,
      py::arg("bars_displace_concrete") = false,
      "Nonlinear service state of a prestressed or reinforced section, each material on its own law (section "
      "arguments and result of prestress_service_state): the concrete on the law dict concrete (kind='sargin' "
      "with fck or fcm, Ecm, eps_c1, eps_cu1; 'concrete_pr'; 'concrete_tension_softening'; 'elastic'), the bonded "
      "tendons on strand from eps_p0 (default: Y1860S7 at its characteristic strengths, gamma_s = 1), the bars on "
      "the bilinear curve Es, fy, ft at eps_ud. The laws are read as given: characteristic or mean strengths for "
      "a service state. A strain past the bound of a law is refused.");
  mod.def(
      "section_prestress_losses",
      [loss_section, loss_dicts](const Rings& rings, const py::list& tendons, double Ecm, double phi, double eps_cs,
                                 const std::vector<double>& relaxation, const Vec3T& sustained,
                                 const py::object& strand) {
        const ec2::RcPolySection s = loss_section(rings, tendons, strand);
        ec2::LossState state;
        state.Ecm = Ecm;
        state.phi = phi;
        state.eps_cs = eps_cs;
        state.relaxation = relaxation;
        state.sustained = Vec3(std::get<0>(sustained), std::get<1>(sustained), std::get<2>(sustained));
        return loss_dicts(engine([&] { return ec2::prestress_losses(s, state); }));
      },
      py::arg("rings"), py::arg("tendons"), py::arg("Ecm"), py::arg("phi"), py::arg("eps_cs"),
      py::arg("relaxation") = std::vector<double>{}, py::arg("sustained") = Vec3T(0.0, 0.0, 0.0),
      py::arg("strand") = py::none(),
      "Time-dependent prestress losses per tendon (EN 1992-1-1 §5.10.6, Expr. (5.46) coupled between the "
      "tendons): creep coefficient phi(t, t0), shrinkage strain eps_cs between t0 and t (negative = shortening), "
      "pure relaxation loss of each tendon [Pa], sustained = (N, My, Mz) besides the prestress. Returns "
      "[dict(stress, elastic_shortening, creep, shrinkage, relaxation, total, eps_p0_after)]; elastic_shortening "
      "= (Ep / Ecm) |sigma_c| under the prestress alone, not in total.");
  mod.def(
      "sequential_tensioning_factor",
      [](int n_tendons) { return engine([&] { return ec2::sequential_tensioning_factor(n_tendons); }); },
      py::arg("n_tendons"),
      "Mean share (n - 1) / (2 n) of the elastic shortening lost by n tendons post-tensioned in succession "
      "(EN 1992-1-1 Expr. (5.44)).");
  mod.def(
      "section_loss_history",
      [service_section, passive, loss_dicts](const Rings& rings, const py::list& tendons,
                                             const std::vector<double>& ages, double fck, double rh, double h0,
                                             double t0, const std::string& cement, double t_shrink0,
                                             bool autogenous, double Ecm, double fpk, int relaxation_class,
                                             double rho1000, const Vec3T& sustained, const py::object& strand,
                                             const Bars& bars, double omega, double Es,
                                             bool bars_displace_concrete) {
        const ec2::RcPolySection s = service_section(rings, tendons, strand, bars, 500e6, 525e6, 0.025, Es);
        ec2::LossLaws laws;
        laws.bars = passive(bars, omega, bars_displace_concrete);
        laws.concrete = concrete_time(fck, rh, h0, cement, t_shrink0, autogenous);
        laws.Ecm = Ecm;
        laws.fpk = fpk;
        laws.relaxation_class = relaxation_class;
        laws.rho1000 = rho1000;
        laws.t0 = t0;
        laws.sustained = Vec3(std::get<0>(sustained), std::get<1>(sustained), std::get<2>(sustained));
        py::list out;
        for (const ec2::LossAtAge& a : engine([&] { return ec2::prestress_loss_history(s, laws, ages); }))
          out.append(py::make_tuple(a.age, loss_dicts(a.losses)));
        return out;
      },
      py::arg("rings"), py::arg("tendons"), py::arg("ages"), py::arg("fck"), py::arg("rh"), py::arg("h0"),
      py::arg("t0") = 28.0, py::arg("cement") = "N", py::arg("t_shrink0") = 3.0, py::arg("autogenous") = true,
      py::arg("Ecm") = 0.0, py::arg("fpk") = 1860e6, py::arg("relaxation_class") = 2, py::arg("rho1000") = 2.5,
      py::arg("sustained") = Vec3T(0.0, 0.0, 0.0), py::arg("strand") = py::none(), py::arg("bars") = Bars{},
      py::arg("omega") = 1.0, py::arg("Es") = 200e9, py::arg("bars_displace_concrete") = false,
      "Losses accumulated from the transfer t0 to each age [days] with the EN 1992-1-1 laws: Annex B creep and "
      "shrinkage (fck, rh [%], h0 [m], cement class, start of drying t_shrink0), §3.3.2 relaxation (fpk, class, "
      "rho1000) on the tendon stress at t0; Ecm = 0 takes the EN value of fck. bars [(y, z, weight)] of area "
      "omega weight: the passive bars counted in the section and in the compatibility (as section_losses). "
      "Returns [(age, [loss dict per tendon])] with the dicts of section_prestress_losses.");

  // ---- time-dependent losses with the passive bars in the compatibility
  mod.def(
      "section_losses",
      [service_section, passive, loss_dicts](const Rings& rings, const py::list& tendons, double Ecm, double phi,
                                             double eps_cs, const std::vector<double>& relaxation,
                                             const Vec3T& sustained, const py::object& strand, const Bars& bars,
                                             double omega, double Es, bool bars_displace_concrete, double ageing,
                                             double relaxation_share) {
        const ec2::RcPolySection s = service_section(rings, tendons, strand, bars, 500e6, 525e6, 0.025, Es);
        ec2::LossState state;
        state.Ecm = Ecm;
        state.phi = phi;
        state.eps_cs = eps_cs;
        state.relaxation = relaxation;
        state.sustained = Vec3(std::get<0>(sustained), std::get<1>(sustained), std::get<2>(sustained));
        state.bars = passive(bars, omega, bars_displace_concrete);
        state.ageing = ageing;
        state.relaxation_share = relaxation_share;
        const ec2::SectionLosses l = engine([&] { return ec2::section_losses(s, state); });
        py::dict d;
        d["tendons"] = loss_dicts(l.tendons);
        d["bars"] = loss_dicts(l.bars);
        d["concrete_change"] = std::make_tuple(l.concrete_change.s0, l.concrete_change.sz, l.concrete_change.sy);
        return d;
      },
      py::arg("rings"), py::arg("tendons"), py::arg("Ecm"), py::arg("phi"), py::arg("eps_cs"),
      py::arg("relaxation") = std::vector<double>{}, py::arg("sustained") = Vec3T(0.0, 0.0, 0.0),
      py::arg("strand") = py::none(), py::arg("bars") = Bars{}, py::arg("omega") = 1.0, py::arg("Es") = 200e9,
      py::arg("bars_displace_concrete") = false, py::arg("ageing") = 0.8, py::arg("relaxation_share") = 0.8,
      "Time-dependent losses of a section with its passive bars (arguments of section_prestress_losses, bars "
      "[(y, z, weight)] of area omega weight and modulus Es): every steel follows the concrete at its level, so "
      "the bars restrain the creep and the shrinkage (age-adjusted effective modulus method; ageing = ageing "
      "coefficient, relaxation_share = share of the pure relaxation kept, both 0.8 in EN 1992-1-1 Expr. (5.46)). "
      "Returns dict(tendons, bars, concrete_change): tendons and bars are lists of the dicts of "
      "section_prestress_losses (for a bar: stress at t0, the tension it loses by cause and eps_p0_after, its "
      "initial strain afterwards — hand it back as the fourth value of the bar, (y, z, weight, eps0), with the "
      "eps_p0_after of the tendons, for the section after its losses); "
      "concrete_change = (s0, sz, sy) is the change of the concrete stress s0 + sz z + sy y between t0 and t "
      "(positive where the compression decreases).");

  // ---- serviceability values of a section (keta/section_sls.h)
  mod.def(
      "stress_limits",
      [](double fck, double fyk, double fpk, double k1, double k2, double k3, double k4, double k5) {
        const ec2::StressLimits l = engine([&] { return ec2::stress_limits(fck, fyk, fpk, {k1, k2, k3, k4, k5}); });
        py::dict d;
        d["concrete_characteristic"] = l.concrete_characteristic;
        d["concrete_quasi_permanent"] = l.concrete_quasi_permanent;
        d["bar"] = l.bar;
        d["bar_imposed"] = l.bar_imposed;
        d["tendon"] = l.tendon;
        return d;
      },
      py::arg("fck"), py::arg("fyk"), py::arg("fpk"), py::arg("k1") = 0.6, py::arg("k2") = 0.45,
      py::arg("k3") = 0.8, py::arg("k4") = 1.0, py::arg("k5") = 0.75,
      "Stress limits of EN 1992-1-1 §7.2 [Pa]: dict(concrete_characteristic = k1 fck, concrete_quasi_permanent = "
      "k2 fck, bar = k3 fyk, bar_imposed = k4 fyk, tendon = k5 fpk), with the recommended factors by default.");
  mod.def(
      "stress_ratios",
      [](const py::dict& state, const py::dict& limits) {
        ec2::PrestressServiceState st;
        st.fibre_min.stress = state["fibre_min"].cast<py::dict>()["stress"].cast<double>();
        st.bar_stress = state["bar_stress"].cast<std::vector<double>>();
        st.tendon_stress = state["tendon_stress"].cast<std::vector<double>>();
        const ec2::StressLimits l{limits["concrete_characteristic"].cast<double>(),
                                  limits["concrete_quasi_permanent"].cast<double>(), limits["bar"].cast<double>(),
                                  limits["bar_imposed"].cast<double>(), limits["tendon"].cast<double>()};
        const ec2::StressRatios r = engine([&] { return ec2::stress_ratios(st, l); });
        py::dict d;
        d["concrete_compression"] = r.concrete_compression;
        d["bar_tension"] = r.bar_tension;
        d["tendon_tension"] = r.tendon_tension;
        d["concrete_characteristic"] = r.concrete_characteristic;
        d["concrete_quasi_permanent"] = r.concrete_quasi_permanent;
        d["bar"] = r.bar;
        d["bar_imposed"] = r.bar_imposed;
        d["tendon"] = r.tendon;
        return d;
      },
      py::arg("state"), py::arg("limits"),
      "The stresses of a service state (a dict returned by prestress_service_state, cracked_service_state, "
      "nonlinear_service_state or service_state) that EN 1992-1-1 §7.2 limits, each over its limit (the dict of "
      "stress_limits): dict(concrete_compression, bar_tension, tendon_tension [Pa], concrete_characteristic, "
      "concrete_quasi_permanent, bar, bar_imposed, tendon [ratios]). Raw values: which ratio applies depends on "
      "the combination the state was computed under.");
  mod.def(
      "service_state",
      [service_section, passive, points_of, state_dict](const Rings& rings, const py::list& tendons, double Ecm,
                                                         double fct_eff, const Vec3T& sustained,
                                                         const py::object& strand, const Points& points,
                                                         const Bars& bars, double omega, double Es,
                                                         bool bars_displace_concrete) {
        const ec2::RcPolySection s = service_section(rings, tendons, strand, bars, 500e6, 525e6, 0.025, Es);
        const ec2::ServiceState r = engine([&] {
          return ec2::service_state(s, Ecm, fct_eff,
                                    Vec3(std::get<0>(sustained), std::get<1>(sustained), std::get<2>(sustained)),
                                    points_of(points), passive(bars, omega, bars_displace_concrete));
        });
        py::dict d = state_dict(r.state);
        d["cracked"] = r.cracked;
        d["flexural_tension"] = r.flexural_tension;
        return d;
      },
      py::arg("rings"), py::arg("tendons"), py::arg("Ecm"), py::arg("fct_eff"),
      py::arg("sustained") = Vec3T(0.0, 0.0, 0.0), py::arg("strand") = py::none(), py::arg("points") = Points{},
      py::arg("bars") = Bars{}, py::arg("omega") = 1.0, py::arg("Es") = 200e9,
      py::arg("bars_displace_concrete") = false,
      "The service state EN 1992-1-1 §7.1(2) selects (arguments of prestress_service_state): the uncracked state "
      "while its largest concrete stress does not exceed fct_eff, the cracked state beyond. Returns the state "
      "dict plus cracked and flexural_tension (the largest concrete stress of the uncracked section). The choice "
      "is made on this action alone.");
  mod.def(
      "section_cracking_moment",
      [service_section, passive](const Rings& rings, const py::list& tendons, double Ecm, double fct, double theta,
                                 const Vec3T& sustained, const py::object& strand, const Bars& bars, double omega,
                                 double Es, bool bars_displace_concrete, bool prestress) {
        ec2::RcPolySection s = service_section(rings, tendons, strand, bars, 500e6, 525e6, 0.025, Es);
        if (!prestress) s = ec2::without_prestress(s);
        const ec2::CrackingMoment c = engine([&] {
          return ec2::cracking_moment(s, Ecm, fct, theta,
                                      Vec3(std::get<0>(sustained), std::get<1>(sustained), std::get<2>(sustained)),
                                      passive(bars, omega, bars_displace_concrete));
        });
        py::dict fibre;
        fibre["y"] = c.fibre.y;
        fibre["z"] = c.fibre.z;
        fibre["stress"] = c.fibre.stress;
        py::dict d;
        d["moment"] = c.moment;
        d["fibre"] = fibre;
        d["modulus"] = c.modulus;
        return d;
      },
      py::arg("rings"), py::arg("tendons"), py::arg("Ecm"), py::arg("fct"), py::arg("theta"),
      py::arg("sustained") = Vec3T(0.0, 0.0, 0.0), py::arg("strand") = py::none(), py::arg("bars") = Bars{},
      py::arg("omega") = 1.0, py::arg("Es") = 200e9, py::arg("bars_displace_concrete") = false,
      py::arg("prestress") = true,
      "Cracking moment of the uncracked section: the moment M_cr of direction theta — (My, Mz) = M_cr (cos "
      "theta, sin theta) — that, added to sustained, brings the first vertex of the rings to fct: (fct - "
      "sigma_cp) W. Returns dict(moment, fibre = dict(y, z, stress = sigma_cp), modulus = W). prestress=False "
      "slackens the tendons (eps_p0 = 0): M_rep of EN 1992-2 §6.1(109) b).");
  mod.def(
      "brittle_failure_steel",
      [](double M_rep, double z, double fyk, double Ap, double fptk) {
        return engine([&] { return ec2::brittle_failure_steel(M_rep, z, fyk, Ap, fptk); });
      },
      py::arg("M_rep"), py::arg("z"), py::arg("fyk"), py::arg("Ap") = 0.0, py::arg("fptk") = 0.0,
      "Minimum reinforcement against brittle failure, EN 1992-2 Expr. (6.101a): M_rep / (z fyk); with tendons "
      "counted (Ap, fptk), Expr. (6.101b): (M_rep / z - Ap min(0.4 fptk, 500 MPa)) / fyk, not below zero [m2].");
  mod.def(
      "ultimate_lever_arm",
      [poly_section](const Rings& rings, double fck, double fcd, const py::list& tendons, const Bars& bars,
                     double omega, double N, double theta, const py::object& strand, double fyd, double ftd,
                     double eps_ud, double Es) {
        const ec2::RcPolySection s = poly_section(rings, bars, fck, fcd, tendons, strand, fyd, ftd, eps_ud, Es);
        return engine([&] { return ec2::ultimate_lever_arm(s, omega, N, theta); });
      },
      py::arg("rings"), py::arg("fck"), py::arg("fcd"), py::arg("tendons"), py::arg("bars") = Bars{},
      py::arg("omega") = 1.0, py::arg("N") = 0.0, py::arg("theta") = 0.0, py::arg("strand") = py::none(),
      py::arg("fyd") = 500e6 / 1.15, py::arg("ftd") = 525e6 / 1.15, py::arg("eps_ud") = 0.025,
      py::arg("Es") = 200e9,
      "Lever arm of the reinforcing steel at the ultimate limit state (arguments of "
      "prestressed_section_resistance): on the ultimate plane carrying N with a moment along theta, the distance "
      "between the resultant of the concrete compression and the resultant of the stretched bars [m].");
  const auto crack_params = [](double fct_eff, double cover, double bar_dia, double spacing, double k1, double k2,
                               double k3, double k4, double kt, double alpha_e, double xi, double tendon_dia,
                               bool xi1_squared) {
    ec2::CrackWidthParams p;
    p.fct_eff = fct_eff;
    p.cover = cover;
    p.bar_dia = bar_dia;
    p.spacing = spacing;
    p.k1 = k1;
    p.k2 = k2;
    p.k3 = k3;
    p.k4 = k4;
    p.kt = kt;
    p.alpha_e = alpha_e;
    p.xi = xi;
    p.tendon_dia = tendon_dia;
    p.xi1_squared = xi1_squared;
    return p;
  };
  const auto crack_dict = [](const ec2::CrackWidth& w) {
    py::dict d;
    d["h"] = w.h;
    d["d"] = w.d;
    d["x"] = w.x;
    d["hc_ef"] = w.hc_ef;
    d["Ac_eff"] = w.Ac_eff;
    d["As"] = w.As;
    d["Ap"] = w.Ap;
    d["xi1"] = w.xi1;
    d["rho_p_eff"] = w.rho_p_eff;
    d["sigma_s"] = w.sigma_s;
    d["d_eps"] = w.d_eps;
    d["s_r_max"] = w.s_r_max;
    d["wide_spacing"] = w.wide_spacing;
    d["wk"] = w.wk;
    return d;
  };
  mod.def(
      "section_crack_width",
      [service_section, passive, crack_params, crack_dict](
          const Rings& rings, const py::list& tendons, double Ecm, const Vec3T& action, double fct_eff,
          const Bars& bars, double omega, double Es, const py::object& strand, bool bars_displace_concrete,
          double cover, double bar_dia, double spacing, double k1, double k2, double k3, double k4, double kt,
          double alpha_e, double xi, double tendon_dia, bool xi1_squared) {
        const ec2::RcPolySection s = service_section(rings, tendons, strand, bars, 500e6, 525e6, 0.025, Es);
        const ec2::CrackWidthParams p =
            crack_params(fct_eff, cover, bar_dia, spacing, k1, k2, k3, k4, kt, alpha_e, xi, tendon_dia, xi1_squared);
        return crack_dict(engine([&] {
          return ec2::crack_width(s, Ecm, Vec3(std::get<0>(action), std::get<1>(action), std::get<2>(action)),
                                  passive(bars, omega, bars_displace_concrete), p);
        }));
      },
      py::arg("rings"), py::arg("tendons"), py::arg("Ecm"), py::arg("action"), py::arg("fct_eff"),
      py::arg("bars") = Bars{}, py::arg("omega") = 1.0, py::arg("Es") = 200e9, py::arg("strand") = py::none(),
      py::arg("bars_displace_concrete") = false, py::arg("cover") = 0.030, py::arg("bar_dia") = 0.016,
      py::arg("spacing") = 0.0, py::arg("k1") = 0.8, py::arg("k2") = 0.5, py::arg("k3") = 3.4,
      py::arg("k4") = 0.425, py::arg("kt") = 0.6, py::arg("alpha_e") = 0.0, py::arg("xi") = 0.0,
      py::arg("tendon_dia") = 0.0, py::arg("xi1_squared") = false,
      "Crack width of a polygonal section under action = (N, My, Mz) (EN 1992-1-1 §7.3.4), from its cracked "
      "state on Ecm: along the strain gradient, hc_ef = min(2.5 (h - d), (h - x) / 3), Ac_eff the concrete within "
      "hc_ef of the most stretched vertex, the bars and bonded tendons of that band in rho_p_eff = (As + xi1 "
      "Ap') / Ac_eff (xi = bond ratio of the tendons, Table 6.2, 0 leaves them out; tendon_dia = their equivalent "
      "diameter; xi1_squared: Expr. (7.10) before the corrigendum AC:2010). alpha_e = 0 takes Es over Ecm. "
      "Returns dict(h, d, x, hc_ef, Ac_eff, As, Ap, xi1, rho_p_eff, sigma_s, d_eps, s_r_max, wide_spacing, wk). A "
      "section without compressed zone is refused: use crack_width_chain with its effective area.");
  mod.def(
      "crack_width_chain",
      [crack_params, crack_dict](double sigma_s, double As, double Ap, double Ac_eff, double h_minus_x, double Es,
                                 double fct_eff, double alpha_e, double cover, double bar_dia, double spacing,
                                 double k1, double k2, double k3, double k4, double kt, double xi,
                                 double tendon_dia, bool xi1_squared) {
        const ec2::CrackWidthParams p =
            crack_params(fct_eff, cover, bar_dia, spacing, k1, k2, k3, k4, kt, alpha_e, xi, tendon_dia, xi1_squared);
        return crack_dict(engine([&] { return ec2::crack_width_of(sigma_s, As, Ap, Ac_eff, h_minus_x, Es, p); }));
      },
      py::arg("sigma_s"), py::arg("As"), py::arg("Ap"), py::arg("Ac_eff"), py::arg("h_minus_x"), py::arg("Es"),
      py::arg("fct_eff"), py::arg("alpha_e"), py::arg("cover") = 0.030, py::arg("bar_dia") = 0.016,
      py::arg("spacing") = 0.0, py::arg("k1") = 0.8, py::arg("k2") = 0.5, py::arg("k3") = 3.4,
      py::arg("k4") = 0.425, py::arg("kt") = 0.6, py::arg("xi") = 0.0, py::arg("tendon_dia") = 0.0,
      py::arg("xi1_squared") = false,
      "EN 1992-1-1 Expr. (7.5) and (7.8) to (7.14) from their terms: steel stress sigma_s of the cracked section, "
      "bars As and bonded tendons Ap within Ac_eff, stretched depth h_minus_x (read by Expr. (7.14)), modulus Es "
      "of the bars. Returns the dict of section_crack_width (h, d, x and hc_ef zero).");

  // ---- composite section: several materials, each on its law over its rings
  const auto composite_of = [law_of](const py::list& parts, const py::list& steels) {
    ec2::CompositeSection s;
    for (const py::handle& h : parts) {
      const py::dict d = py::reinterpret_borrow<py::dict>(h);
      ec2::MaterialPart part;
      for (const auto& r : d["rings"].cast<Rings>()) {
        std::vector<Vec2> ring;
        for (const auto& [y, z] : r) ring.emplace_back(y, z);
        part.rings.push_back(ring);
      }
      part.law = law_of(d["law"].cast<py::dict>(), "composite section");
      s.parts.push_back(part);
    }
    for (const py::handle& h : steels) {
      const py::dict d = py::reinterpret_borrow<py::dict>(h);
      ec2::PointSteel t;
      t.y = d["y"].cast<double>();
      t.z = d["z"].cast<double>();
      t.A = d["A"].cast<double>();
      t.law = law_of(d["law"].cast<py::dict>(), "composite section");
      if (d.contains("eps0")) t.eps0 = d["eps0"].cast<double>();
      s.steels.push_back(t);
    }
    return s;
  };
  mod.def(
      "composite_section_resultants",
      [composite_of](const py::list& parts, const Vec3T& plane, const py::list& steels) {
        const ec2::CompositeSection s = composite_of(parts, steels);
        const Vec3 r = engine([&] {
          return ec2::composite_resultants(s, {std::get<0>(plane), std::get<1>(plane), std::get<2>(plane)});
        });
        return Vec3T(r(0), r(1), r(2));
      },
      py::arg("parts"), py::arg("plane"), py::arg("steels") = py::list(),
      "(N, My, Mz) of a composite section under the plane (e0, ky, kz), strain e0 + ky z - kz y: parts = "
      "[dict(rings, law)], each material on its own law dict over its own rings (outer CCW, voids CW; a part "
      "inside another is a void of the outer one), steels = [dict(y, z, A, law, eps0=0)] point steels. No mesh: "
      "contour integrals of each law. Laws of a part: 'concrete_pr', 'sargin', 'concrete_tension_softening', "
      "'elastic', 'steel', 'steel_hardening'; of a point steel: 'elastic', 'steel', 'steel_hardening', 'strand'.");
  mod.def(
      "composite_section_state",
      [composite_of](const py::list& parts, const Vec3T& action, const py::list& steels) {
        const ec2::CompositeSection s = composite_of(parts, steels);
        const ec2::CompositeState st = engine([&] {
          return ec2::composite_state(s, Vec3(std::get<0>(action), std::get<1>(action), std::get<2>(action)));
        });
        std::vector<Vec3T> carried;
        for (const Vec3& r : st.part_resultants) carried.emplace_back(r(0), r(1), r(2));
        py::dict d;
        d["plane"] = std::make_tuple(st.plane.e0, st.plane.ky, st.plane.kz);
        d["part_resultants"] = carried;
        d["steel_stress"] = st.steel_stress;
        return d;
      },
      py::arg("parts"), py::arg("action"), py::arg("steels") = py::list(),
      "State of a composite section (arguments of composite_section_resultants) under action = (N, My, Mz): the "
      "plane whose resultants are the action. Returns dict(plane = (e0, ky, kz), part_resultants = [(N, My, Mz) "
      "carried by each part], steel_stress). A state straining a part or a steel past the bound of its law is "
      "refused.");
}
