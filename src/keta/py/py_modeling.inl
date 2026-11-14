// keta Python bindings — modelling and result completions: end
// eccentricities, polygonal and thin-walled sections, nodal and bar masses,
// imposed strains, shell temperatures, global activation, tendons in a static
// case, EN 1990 combination generation, envelopes with concomitant states,
// EN 1998-1 spectra. Included by keta_py.cpp inside its anonymous namespace.
//
//   m.set_eccentricity(b, origin=(0, 0, 0.3), end=(0, 0, 0.5))
//   s = m.add_section_polygon([(-0.5, -0.3), (0.5, -0.3), (0.5, 0.3), (-0.5, 0.3)])
//   m.add_prestrain(b, strain=-2e-4, ky=1e-4, case=3)
//   lo, hi = m.analyze().envelope(cases=[1, 2, 3])

Vec3 vec3_of(const Vec3T& a) { return Vec3(std::get<0>(a), std::get<1>(a), std::get<2>(a)); }

std::vector<Vec2> ring_of(const std::vector<std::pair<double, double>>& pts, const char* who) {
  if (pts.size() < 3) throw py::value_error(std::string(who) + ": a ring needs at least 3 points (y, z)");
  std::vector<Vec2> r;
  for (const auto& [y, z] : pts) {
    if (!std::isfinite(y) || !std::isfinite(z)) throw py::value_error(std::string(who) + ": non-finite coordinate");
    r.emplace_back(y, z);
  }
  return r;
}

py::dict exact_dict(const ExactSection& e) {
  py::dict d;
  d["A"] = e.A;
  d["yc"] = e.yc;
  d["zc"] = e.zc;
  d["Iy"] = e.Iy;
  d["Iz"] = e.Iz;
  d["Iyz"] = e.Iyz;
  d["I1"] = e.I1;
  d["I2"] = e.I2;
  d["phi"] = e.phi;
  d["J"] = e.J;
  d["Iw"] = e.Iw;
  d["ys"] = e.y_sc;
  d["zs"] = e.z_sc;
  d["Avy"] = e.Avy;
  d["Avz"] = e.Avz;
  d["ymax"] = e.ymax;
  d["zmax"] = e.zmax;
  d["EA_ref"] = e.EA_ref;
  d["E_ref"] = e.E_ref;
  return d;
}

ExactInput exact_input(const std::vector<std::pair<double, double>>& outer,
                       const std::vector<std::vector<std::pair<double, double>>>& holes, double nu, double mesh_size) {
  ExactInput in;
  in.outline = ring_of(outer, "section");
  for (const auto& h : holes) in.holes.push_back(ring_of(h, "section hole"));
  ExactMaterial mat;
  mat.nu = nu;
  in.materials = {mat};
  in.target_h = mesh_size;
  return in;
}

// The materials [dict(E, nu=0.3)] and the regions [dict(outline, material)] of a
// composite section: material 0 is the reference and fills what no region covers.
void exact_composite(ExactInput& in, const py::list& materials, const py::list& regions) {
  if (py::len(materials) == 0) {
    if (py::len(regions) != 0) throw py::value_error("section: regions need materials");
    return;
  }
  in.materials.clear();
  for (const py::handle& h : materials) {
    const py::dict d = py::reinterpret_borrow<py::dict>(h);
    if (!d.contains("E")) throw py::value_error("section material: dict(E, nu=0.3)");
    ExactMaterial m;
    m.E = d["E"].cast<double>();
    if (d.contains("nu")) m.nu = d["nu"].cast<double>();
    in.materials.push_back(m);
  }
  for (const py::handle& h : regions) {
    const py::dict d = py::reinterpret_borrow<py::dict>(h);
    if (!d.contains("outline") || !d.contains("material")) throw py::value_error("section region: dict(outline, material)");
    ExactRegion r;
    r.outline = ring_of(d["outline"].cast<std::vector<std::pair<double, double>>>(), "section region");
    r.material = d["material"].cast<int>();
    if (r.material < 0 || r.material >= static_cast<int>(in.materials.size()))
      throw py::value_error("section region: material is an index of materials");
    in.regions.push_back(r);
  }
}

// A result-like bound (envelope bound or concomitant state): every accessor of
// CaseResult. `c` may alias into a StaticResults it keeps alive.
struct PyBound {
  ModelPtr model;
  std::shared_ptr<const CaseResult> c;
  std::string label;
};

std::vector<const CaseResult*> sources_of(const PyResults& r, const std::optional<std::vector<Id>>& cases,
                                          const std::optional<std::vector<Id>>& combos, std::vector<std::string>* names) {
  std::vector<const CaseResult*> v;
  const bool all = !cases && !combos;
  if (all || cases) {
    for (const CaseResult& c : r.r->cases)
      if (all || std::find(cases->begin(), cases->end(), c.id) != cases->end()) {
        v.push_back(&c);
        if (names) names->push_back("case " + std::to_string(c.id));
      }
    if (cases)
      for (Id id : *cases)
        if (std::none_of(r.r->cases.begin(), r.r->cases.end(), [&](const CaseResult& c) { return c.id == id; }))
          throw py::key_error("no load case " + std::to_string(id));
  }
  if (all || combos) {
    for (const CaseResult& c : r.r->combos)
      if (all || std::find(combos->begin(), combos->end(), c.id) != combos->end()) {
        v.push_back(&c);
        if (names) names->push_back("combination " + std::to_string(c.id));
      }
    if (combos)
      for (Id id : *combos)
        if (std::none_of(r.r->combos.begin(), r.r->combos.end(), [&](const CaseResult& c) { return c.id == id; }))
          throw py::key_error("no combination " + std::to_string(id));
  }
  if (v.empty()) throw py::value_error("envelope: no source result");
  return v;
}

VariableAction variable_of(const py::handle& h) {
  if (py::isinstance<py::tuple>(h) || py::isinstance<py::list>(h)) {
    const auto [id, kind] = h.cast<std::pair<Id, std::string>>();
    if (kind == "traffic") return traffic_gr1a(id);
    if (kind == "wind") return wind_action(id);
    if (kind == "thermal") return thermal_action(id);
    if (kind == "snow") return snow_action(id);
    throw py::value_error("variable action kind '" + kind + "': 'traffic', 'wind', 'thermal' or 'snow'");
  }
  const py::dict d = py::reinterpret_borrow<py::dict>(h);
  VariableAction a;
  a.case_id = d["case"].cast<Id>();
  a.name = d.contains("name") ? d["name"].cast<std::string>() : "Q" + std::to_string(a.case_id);
  if (d.contains("psi0")) a.psi0 = d["psi0"].cast<double>();
  if (d.contains("psi1")) a.psi1 = d["psi1"].cast<double>();
  if (d.contains("psi2")) a.psi2 = d["psi2"].cast<double>();
  if (d.contains("gamma_q")) a.gamma_q = d["gamma_q"].cast<double>();
  return a;
}

void bind_modeling(py::module_& mod) {
  auto M = py::reinterpret_borrow<py::class_<PyModel>>(mod.attr("Model"));

  // ---------------- geometry ----------------
  M.def(
      "set_eccentricity",
      [](PyModel& self, Id bar, Vec3T origin, Vec3T end) {
        Bar& b = self.bar(bar);
        b.ecc = {vec3_of(origin), vec3_of(end)};
      },
      py::arg("bar"), py::arg("origin") = Vec3T{0, 0, 0}, py::arg("end") = Vec3T{0, 0, 0},
      "Rigid-arm offsets of the two ends of a bar, global frame [m]: the bar runs from node_i + origin to "
      "node_j + end.");

  M.def(
      "add_section_polygon",
      [](PyModel& self, const std::vector<std::pair<double, double>>& outer,
         const std::vector<std::vector<std::pair<double, double>>>& holes, std::optional<std::string> name, bool exact,
         double nu, double mesh_size) {
        Section s;
        s.name = name.value_or("Polygon");
        ShapePolygon sh;
        for (const auto& [y, z] : outer) sh.outer.emplace_back(y, z, 0.0);
        for (const auto& h : holes) {
          std::vector<Vec3> r;
          for (const auto& [y, z] : h) r.emplace_back(y, z, 0.0);
          sh.holes.push_back(r);
        }
        (void)ring_of(outer, "add_section_polygon");
        s.shape = sh;
        engine([&] {
          s.finalize();
          if (exact) attach_exact(s, exact_input(outer, holes, nu, mesh_size));
          return 0;
        });
        s.id = self.next_section++;
        self.m.sections.push_back(std::move(s));
        self.m.invalidate_index();
        return self.m.sections.back().id;
      },
      py::arg("outer"), py::arg("holes") = std::vector<std::vector<std::pair<double, double>>>{},
      py::arg("name") = py::none(), py::arg("exact") = true, py::arg("nu") = 0.2, py::arg("mesh_size") = 0.0,
      "Section of any polygonal shape: outer = [(y, z), ...] [m], holes = list of rings. Area and second moments by "
      "Green's theorem (the axes must be principal); with exact=True the torsion constant and the shear areas come "
      "from a quadratic finite-element solution of the Saint-Venant warping problem (Poisson ratio nu, element "
      "size mesh_size, 0 = automatic).");

  mod.def(
      "section_polygon",
      [](const std::vector<std::pair<double, double>>& outer,
         const std::vector<std::vector<std::pair<double, double>>>& holes, double nu, double mesh_size,
         const py::list& materials, const py::list& regions) {
        ExactInput in = exact_input(outer, holes, nu, mesh_size);
        exact_composite(in, materials, regions);
        return exact_dict(engine([&] { return analyze_exact(in); }));
      },
      py::arg("outer"), py::arg("holes") = std::vector<std::vector<std::pair<double, double>>>{}, py::arg("nu") = 0.2,
      py::arg("mesh_size") = 0.0, py::arg("materials") = py::list(), py::arg("regions") = py::list(),
      "Properties of a polygonal section with openings (Tri6 finite elements of the warping problem): dict(A, yc, "
      "zc, Iy, Iz, Iyz (centroidal), I1, I2, phi (principal), J, Iw, ys, zs (shear centre, from the centroid), Avy, "
      "Avz, ymax, zmax, EA_ref, E_ref) [m]. Composite section: materials = [dict(E, nu=0.3)] and regions = "
      "[dict(outline, material)] — the elements inside a region take its material (an index of materials), the "
      "others material 0; the properties are those of the section transformed to E_ref = the modulus of material "
      "0 (A stays the geometric area, EA_ref is the transformed one), and nu is then read from the materials.");

  mod.def(
      "section_stress",
      [](const std::vector<std::pair<double, double>>& outer, const py::dict& actions,
         const std::vector<std::pair<double, double>>& points,
         const std::vector<std::vector<std::pair<double, double>>>& holes, double nu, double mesh_size,
         const py::list& materials, const py::list& regions) {
        ExactInput in = exact_input(outer, holes, nu, mesh_size);
        exact_composite(in, materials, regions);
        SectionActions a;
        for (const auto& [key, value] : actions) {
          const std::string name = py::str(key);
          const double v = value.cast<double>();
          if (name == "N") a.N = v;
          else if (name == "Vy") a.Vy = v;
          else if (name == "Vz") a.Vz = v;
          else if (name == "T") a.T = v;
          else if (name == "My") a.My = v;
          else if (name == "Mz") a.Mz = v;
          else if (name == "B") a.B = v;
          else throw py::value_error("section_stress: actions N, Vy, Vz, T, My, Mz, B");
        }
        const ExactSection e = engine([&] { return analyze_exact(in); });
        py::list at;
        for (const auto& [y, z] : points) {
          const StressState s = engine([&] { return e.stress_at(y - e.yc, z - e.zc, a); });
          py::dict d;
          d["sigma"] = s.sig;
          d["tau_y"] = s.tau_y;
          d["tau_z"] = s.tau_z;
          d["tau"] = s.tau;
          d["von_mises"] = s.vm;
          at.append(d);
        }
        const SigmaRange r = engine([&] { return e.sigma_extremes(a); });
        py::dict out;
        out["points"] = at;
        out["sigma_min"] = Vec3T(r.min.value, r.min.y + e.yc, r.min.z + e.zc);
        out["sigma_max"] = Vec3T(r.max.value, r.max.y + e.yc, r.max.z + e.zc);
        out["von_mises_max"] = engine([&] { return e.von_mises_max(a); });
        return out;
      },
      py::arg("outer"), py::arg("actions"), py::arg("points") = std::vector<std::pair<double, double>>{},
      py::arg("holes") = std::vector<std::vector<std::pair<double, double>>>{}, py::arg("nu") = 0.2,
      py::arg("mesh_size") = 0.0, py::arg("materials") = py::list(), py::arg("regions") = py::list(),
      "Stresses of a polygonal section (arguments of section_polygon) under the actions dict(N, Vy, Vz, T, My, Mz, "
      "B) taken at its elastic centroid and shear centre [N, N.m, N.m^2]: sigma = N / A + My z / Iy - Mz y / Iz "
      "(times E / E_ref in a composite section) plus the warping stress of B, shear stresses of Vy, Vz and T from "
      "the Saint-Venant fields. points = [(y, z)] in the frame of the outline. Returns dict(points = [dict(sigma, "
      "tau_y, tau_z, tau, von_mises)], sigma_min, sigma_max = (value, y, z) of the normal stress of N, My, Mz, "
      "von_mises_max over the mesh nodes) [Pa].");

  mod.def(
      "import_dxf",
      [](const std::string& text, double scale, double arc_sagitta, bool is_path) {
        DxfOptions o;
        o.scale = scale;
        o.arc_sagitta = arc_sagitta;
        const DxfSection s = engine([&] { return is_path ? load_dxf(text, o) : import_dxf(text, o); });
        const auto ring = [](const std::vector<Vec2>& r) {
          std::vector<std::pair<double, double>> out;
          for (const Vec2& p : r) out.emplace_back(p.x(), p.y());
          return out;
        };
        std::vector<std::vector<std::pair<double, double>>> holes;
        for (const auto& h : s.holes) holes.push_back(ring(h));
        py::dict d;
        d["outer"] = ring(s.outline);
        d["holes"] = holes;
        return d;
      },
      py::arg("text"), py::arg("scale") = 1.0, py::arg("arc_sagitta") = 0.0, py::arg("is_path") = false,
      "Outline of a section drawn in a DXF text (ASCII; closed LWPOLYLINE / POLYLINE with their arcs, CIRCLE), or "
      "in the file `text` when is_path: the largest closed contour is the outer boundary, the others its "
      "openings. scale converts the drawing unit to metres (1e-3 for millimetres); arc_sagitta is the largest "
      "chord-to-arc distance of the polygonised arcs in drawing units (0: one thousandth of each radius). Returns "
      "dict(outer, holes), counter-clockwise, ready for section_polygon and add_section_polygon.");

  mod.def(
      "thinwall_section",
      [](const std::vector<std::pair<double, double>>& nodes, const py::list& walls) {
        ThinWall w;
        for (const auto& [y, z] : nodes) w.nodes.emplace_back(y, z);
        for (const py::handle& h : walls) {
          const py::tuple wall = py::reinterpret_borrow<py::tuple>(h);
          if (wall.size() != 3 && wall.size() != 4)
            throw py::value_error("thinwall_section: a wall is (node_a, node_b, t) or (node_a, node_b, t_a, t_b)");
          WallSegment seg{wall[0].cast<int>(), wall[1].cast<int>(), wall[2].cast<double>()};
          if (wall.size() == 4) seg.t_b = wall[3].cast<double>();
          w.segs.push_back(seg);
        }
        const ThinWallProps p = engine([&] {
          w.validate();
          return thinwall_props(w);
        });
        py::dict d;
        d["A"] = p.A;
        d["yc"] = p.yc;
        d["zc"] = p.zc;
        d["Iy"] = p.Iy;
        d["Iz"] = p.Iz;
        d["Iyz"] = p.Iyz;
        d["I1"] = p.I1;
        d["I2"] = p.I2;
        d["phi"] = p.phi;
        d["J"] = p.J;
        d["J_open"] = p.J_open;
        d["J_closed"] = p.J_closed;
        d["cells"] = p.n_cells;
        d["ys"] = p.y_sc;
        d["zs"] = p.z_sc;
        d["Iw"] = p.Iw;
        d["Avy"] = p.Avy;
        d["Avz"] = p.Avz;
        d["Wely"] = p.Wely;
        d["Welz"] = p.Welz;
        return d;
      },
      py::arg("nodes"), py::arg("walls"),
      "Thin-walled section from its midline: nodes = [(y, z)], walls = [(node_a, node_b, thickness)] (indices from "
      "0), or (node_a, node_b, t_a, t_b) for a wall whose thickness varies linearly from node a to node b. Open "
      "parts by Vlasov, closed cells by Bredt (any number of cells): dict(A, yc, zc, Iy, Iz, Iyz, I1, I2, phi, J, "
      "J_open, J_closed, cells, ys, zs, Iw, Avy, Avz, Wely, Welz) [m].");

  // the wall graph as thinwall_section takes it: nodes [(y, z)], walls [(node_a, node_b, t_a, t_b)]
  const auto wall_graph_dict = [](const ThinWall& w) {
    py::dict d;
    py::list nodes, walls;
    for (const Vec2& p : w.nodes) nodes.append(py::make_tuple(p.x(), p.y()));
    for (const WallSegment& s : w.segs) walls.append(py::make_tuple(s.a, s.b, s.t, s.thickness_at_b()));
    d["nodes"] = nodes;
    d["walls"] = walls;
    return d;
  };
  const auto contour_points = [](const std::vector<std::pair<double, double>>& contour) {
    std::vector<Vec2> points;
    for (const auto& [y, z] : contour) points.emplace_back(y, z);
    return points;
  };
  mod.def(
      "thinwall_from_quads",
      [wall_graph_dict, contour_points](const std::vector<std::pair<double, double>>& contour,
                                        const std::vector<std::array<int, 4>>& walls, double snap) {
        return wall_graph_dict(engine([&] { return thinwall_from_quads(contour_points(contour), walls, snap); }));
      },
      py::arg("contour"), py::arg("walls"), py::arg("snap") = 0.0,
      "Wall graph of a section from its contour points [(y, z)]: each wall is four point indices (p, q, r, s) in "
      "order around it, p->q and s->r its two faces. The midline joins the middles of the end caps p-s and q-r, "
      "the thickness of an end is its cap projected on the normal of the midline; cap middles closer than snap "
      "(0: 1e-6 of the extent) merge. Returns dict(nodes, walls) for thinwall_section, walls as (node_a, node_b, "
      "t_a, t_b).");
  mod.def(
      "thinwall_from_nodes",
      [wall_graph_dict, contour_points](const std::vector<std::pair<double, double>>& contour,
                                        const std::vector<std::vector<int>>& nodes,
                                        const std::vector<std::array<int, 2>>& branches) {
        return wall_graph_dict(engine([&] { return thinwall_from_nodes(contour_points(contour), nodes, branches); }));
      },
      py::arg("contour"), py::arg("nodes"), py::arg("branches"),
      "Wall graph of a section from its contour points [(y, z)], the contour points gathered at each junction "
      "(nodes: lists of at least two point indices) and the branches (node_a, node_b). The spread of a node's "
      "points across a branch is its thickness there; each node lands on the intersection of the midlines of its "
      "branches. Returns dict(nodes, walls) for thinwall_section, walls as (node_a, node_b, t_a, t_b).");

  mod.def(
      "resal_web_shear",
      [](double A, double Iy, double Iz, const std::vector<std::tuple<double, double, double, double, double>>& chords,
         double N, double Vy, double Vz, double My, double Mz) {
        ResalSection s;
        s.A = A;
        s.Iy = Iy;
        s.Iz = Iz;
        for (const auto& [area, y, z, dy_dx, dz_dx] : chords) s.chords.push_back({area, y, z, dy_dx, dz_dx});
        const ResalFactors k = engine([&] { return resal_factors(s); });
        const WebShear w = engine([&] { return resal_web_shear(k, N, Vy, Vz, My, Mz); });
        py::dict d;
        d["Vy"] = w.Vy;
        d["Vz"] = w.Vz;
        d["relief_y"] = w.relief_y;
        d["relief_z"] = w.relief_z;
        d["factors_y"] = std::make_tuple(k.y(0), k.y(1), k.y(2));
        d["factors_z"] = std::make_tuple(k.z(0), k.z(1), k.z(2));
        return d;
      },
      py::arg("A"), py::arg("Iy"), py::arg("Iz"), py::arg("chords"), py::arg("N") = 0.0, py::arg("Vy") = 0.0,
      py::arg("Vz") = 0.0, py::arg("My") = 0.0, py::arg("Mz") = 0.0,
      "Résal effect of a variable-depth section (A, Iy, Iz about its principal centroidal axes): chords = [(area, "
      "y, z, dy_dx, dz_dx)], centroid from the centroid of the section and slopes of its line along the bar. Each "
      "chord carries sigma_f area along its line; returns dict(Vy, Vz: shear left to the web, relief_y, relief_z: "
      "shear carried by the chords, factors_y, factors_z: relief per unit of (N, My, Mz)).");

  // ---------------- masses ----------------
  M.def(
      "add_nodal_mass",
      [](PyModel& self, Id node, std::optional<std::vector<double>> diag, std::optional<std::vector<double>> matrix,
         bool add) {
        (void)self.m.node(node);
        NodalMass nm;
        nm.node = node;
        nm.add = add;
        if (diag.has_value() == matrix.has_value())
          throw py::value_error("add_nodal_mass: give diag= (6 values) or matrix= (36 values or 21 lower-triangle)");
        if (diag) {
          if (diag->size() != 6) throw py::value_error("add_nodal_mass: diag needs 6 values (m, m, m, Ix, Iy, Iz)");
          for (int i = 0; i < 6; ++i) nm.m(i, i) = (*diag)[static_cast<std::size_t>(i)];
        } else if (matrix->size() == 36) {
          for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j) nm.m(i, j) = (*matrix)[static_cast<std::size_t>(6 * i + j)];
          if (!nm.m.isApprox(nm.m.transpose(), 1e-12)) throw py::value_error("add_nodal_mass: matrix not symmetric");
        } else if (matrix->size() == 21) {
          std::size_t k = 0;
          for (int i = 0; i < 6; ++i)
            for (int j = 0; j <= i; ++j) nm.m(i, j) = nm.m(j, i) = (*matrix)[k++];
        } else {
          throw py::value_error("add_nodal_mass: matrix needs 36 values or 21 (lower triangle by rows)");
        }
        self.m.nodal_masses.push_back(nm);
      },
      py::arg("node"), py::arg("diag") = py::none(), py::arg("matrix") = py::none(), py::arg("add") = false,
      "Nodal mass matrix in the global frame (dx dy dz rx ry rz) [kg, kg m^2]: diag = 6 terms, or matrix = 36 terms "
      "(rows) or 21 (lower triangle by rows). add=False: it replaces the matrix given earlier at that node; "
      "add=True: it is added to it. The automatic masses of the bars are not affected.");

  M.def(
      "add_bar_mass",
      [](PyModel& self, Id bar, std::string kind, double m1, double m2, double m3, std::optional<Vec3T> ecc_origin,
         std::optional<Vec3T> ecc_end) {
        (void)self.bar(bar);
        BarMass b;
        b.bar = bar;
        for (char& c : kind) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (kind == "uniform") b.kind = BarMassKind::Uni;
        else if (kind == "linear") b.kind = BarMassKind::Lin;
        else if (kind == "parabolic") b.kind = BarMassKind::Parabol;
        else if (kind == "none") b.kind = BarMassKind::Nul;
        else if (kind == "auto") b.kind = BarMassKind::Auto;
        else throw py::value_error("add_bar_mass: kind 'uniform', 'linear', 'parabolic', 'none' or 'auto'");
        if (!(m1 >= 0 && m2 >= 0 && m3 >= 0)) throw py::value_error("add_bar_mass: line masses must be >= 0");
        b.m1 = m1;
        b.m2 = kind == "uniform" ? m1 : m2;
        b.m3 = m3;
        if (ecc_origin) b.ecc_or = vec3_of(*ecc_origin);
        if (ecc_end) b.ecc_ex = vec3_of(*ecc_end);
        self.m.bar_masses.push_back(b);
      },
      py::arg("bar"), py::arg("kind") = "uniform", py::arg("m1") = 0.0, py::arg("m2") = 0.0, py::arg("m3") = 0.0,
      py::arg("ecc_origin") = py::none(), py::arg("ecc_end") = py::none(),
      "Line mass ADDED to a bar [kg/m] for the dynamic analyses: 'uniform' (m1), 'linear' (m1 at the origin, m2 at "
      "the end), 'parabolic' (m1, m2 at mid-length, m3 at the end), offset from the bar axis by ecc_origin / ecc_end "
      "(global frame [m]). 'none' removes the bar's automatic rho A mass; 'auto' keeps it (with "
      "set_bar_mass_selective).");

  M.def(
      "set_bar_mass_selective",
      [](PyModel& self, bool on) { self.m.bar_mass_selective = on; }, py::arg("on") = true,
      "on=True: only the bars given an add_bar_mass(kind='auto') entry keep their automatic rho A mass.");

  // ---------------- loads ----------------
  M.def(
      "add_prestrain",
      [](PyModel& self, Id bar, double strain, double ky, double kz, Id cs, double kx) {
        (void)self.bar(bar);
        Load l = Load::prestrain(bar, strain, ky, kz);
        l.kappa_x = kx;
        self.case_mut(cs).loads.push_back(l);
      },
      py::arg("bar"), py::arg("strain") = 0.0, py::arg("ky") = 0.0, py::arg("kz") = 0.0, py::arg("case") = 1,
      py::arg("kx") = 0.0,
      "Uniform imposed strain on a bar: axial strain (> 0 lengthening) and curvatures ky, kz [1/m] about the local "
      "axes, the opposites of the rotation gradients: ky = -d(theta_y)/dx = d2w/dx2 (> 0 lengthens the fibre -z, as "
      "a thermal gradient gz = -ky / alpha) and kz = -d(theta_z)/dx = -d2v/dx2 (> 0 lengthens the fibre +y, as "
      "gy = kz / alpha), and the twist rate kx = d(theta_x)/dx [rad/m] (the free twist of the bar).");
  M.def(
      "add_prestrain_linear",
      [](PyModel& self, Id bar, Vec3T start, Vec3T end, Id cs, std::pair<double, double> kx) {
        (void)self.bar(bar);
        Load l = Load::prestrain_linear(bar, std::get<0>(start), std::get<1>(start), std::get<2>(start),
                                        std::get<0>(end), std::get<1>(end), std::get<2>(end));
        l.kappa_x = kx.first;
        l.kappa_x2 = kx.second;
        self.case_mut(cs).loads.push_back(l);
      },
      py::arg("bar"), py::arg("start"), py::arg("end"), py::arg("case") = 1,
      py::arg("kx") = std::pair<double, double>(0.0, 0.0),
      "Imposed strain varying linearly along the bar: start = (strain, ky, kz) at the origin, end = the same at the "
      "end; kx = (start, end) twist rate d(theta_x)/dx [rad/m].");
  M.def(
      "add_prestrain_parabolic",
      [](PyModel& self, Id bar, Vec3T start, Vec3T mid, Vec3T end, Id cs, std::tuple<double, double, double> kx) {
        (void)self.bar(bar);
        auto v = [](const Vec3T& t) { return Vec3(std::get<0>(t), std::get<1>(t), std::get<2>(t)); };
        Load l = Load::prestrain_parabolic(bar, v(start), v(mid), v(end));
        l.kappa_x = std::get<0>(kx);
        l.kappa_x_m = std::get<1>(kx);
        l.kappa_x2 = std::get<2>(kx);
        self.case_mut(cs).loads.push_back(l);
      },
      py::arg("bar"), py::arg("start"), py::arg("mid"), py::arg("end"), py::arg("case") = 1,
      py::arg("kx") = std::tuple<double, double, double>(0.0, 0.0, 0.0),
      "Imposed strain varying as a parabola along the bar: start, mid and end = (strain, ky, kz) at the origin, at "
      "mid-bar and at the end (the parabola through them; the creep curvature of a bar under a distributed load); "
      "kx = (start, mid, end) twist rate d(theta_x)/dx [rad/m]. Signs as add_prestrain.");
  M.def(
      "add_prestrain_point",
      [](PyModel& self, Id bar, double at, double dx, double ry, double rz, double dy, double dz, double rx,
         bool absolute, Id cs) {
        (void)self.bar(bar);
        const double L = bar_length(self.m, bar);
        const double a = absolute ? at : at * L;
        if (!(a >= 0 && a <= L * (1 + 1e-12))) throw py::value_error("add_prestrain_point: abscissa outside the bar");
        Load l = Load::prestrain_point(bar, a, dx, ry, rz);
        l.dvy = dy;
        l.dvz = dz;
        l.dthx = rx;
        self.case_mut(cs).loads.push_back(l);
      },
      py::arg("bar"), py::arg("at"), py::arg("dx") = 0.0, py::arg("ry") = 0.0, py::arg("rz") = 0.0,
      py::arg("dy") = 0.0, py::arg("dz") = 0.0, py::arg("rx") = 0.0, py::arg("absolute") = false, py::arg("case") = 1,
      "Concentrated dislocation in a bar at `at` (relative, or [m] with absolute=True), local frame: gaps dx, dy, dz "
      "[m] and kinks rx, ry, rz [rad] between the two sides of the section.");
  M.def(
      "add_surface_thermal",
      [](PyModel& self, Id slab, double dt, double grad, Id cs) {
        self.case_mut(cs).loads.push_back(Load::surface_thermal(slab, dt, grad));
      },
      py::arg("slab"), py::arg("dt"), py::arg("grad") = 0.0, py::arg("case") = 1,
      "Temperature of the shells of a slab: uniform change dt [K] and gradient grad = dT/dz through the thickness "
      "[K/m] (along the facet normal). The material needs alpha.");

  // ---------------- activation ----------------
  M.def(
      "set_active",
      [](PyModel& self, std::optional<std::vector<Id>> bars, std::optional<std::vector<Id>> supports) {
        if (bars) {
          for (Id b : *bars) (void)self.bar(b);
          self.active_bars = std::set<Id>(bars->begin(), bars->end());
        } else {
          self.active_bars.reset();
        }
        if (supports) {
          for (Id n : *supports) (void)self.m.node(n);
          self.active_supports = std::set<Id>(supports->begin(), supports->end());
        } else {
          self.active_supports.reset();
        }
      },
      py::arg("bars") = py::none(), py::arg("supports") = py::none(),
      "Global activation for the analyses that follow: bars = the active bars, supports = the nodes whose supports "
      "are active (None = all). Loads and added masses on inactive bars are dropped; a node left without any element "
      "is held still. Phased analyses ignore it (they have their own activation).");

  // ---------------- tendons in a static case ----------------
  M.def(
      "add_tendon_loads",
      [](PyModel& self, Id tendon, Id cs, double fraction) {
        for (Tendon& t : self.m.tendons)
          if (t.id == tendon) {
            ModelPtr model = self.solver_model();
            const std::vector<Load> loads = engine([&] { return tendon_loads(*model, t, fraction); });
            LoadCase& c = self.case_mut(cs);
            c.loads.insert(c.loads.end(), loads.begin(), loads.end());
            t.case_id = cs;
            return static_cast<int>(loads.size());
          }
        throw py::key_error("no tendon " + std::to_string(tendon));
      },
      py::arg("tendon"), py::arg("case") = 1, py::arg("fraction") = 1.0,
      "Puts the equivalent loads of a tendon (add_tendon) in a static load case: its force after friction and "
      "draw-in times `fraction`, without the time-dependent losses. Returns the number of loads added. Not for "
      "phased analyses, which tension the tendons themselves.");

  // ---------------- incremental launching (keta/launching.h) ----------------
  const auto launch_path = [](const std::string& kind, const Vec3T& origin, double heading, double radius,
                              double grade) {
    LaunchPath p;
    if (kind == "straight") p.kind = LaunchPath::Kind::Straight;
    else if (kind == "circular") p.kind = LaunchPath::Kind::Circular;
    else if (kind == "helical") p.kind = LaunchPath::Kind::Helical;
    else throw py::value_error("launch path: kind 'straight', 'circular' or 'helical'");
    p.origin = Vec3(std::get<0>(origin), std::get<1>(origin), std::get<2>(origin));
    p.heading = heading;
    p.radius = radius;
    p.grade = grade;
    return p;
  };
  mod.def(
      "launch_path_point",
      [launch_path](double s, const std::string& kind, const Vec3T& origin, double heading, double radius,
                    double grade) {
        const LaunchPoint p = engine([&] { return launch_path_point(launch_path(kind, origin, heading, radius, grade), s); });
        return std::make_pair(Vec3T(p.x.x(), p.x.y(), p.x.z()), Vec3T(p.t.x(), p.t.y(), p.t.z()));
      },
      py::arg("s"), py::arg("kind") = "straight", py::arg("origin") = Vec3T(0.0, 0.0, 0.0), py::arg("heading") = 0.0,
      py::arg("radius") = 0.0, py::arg("grade") = 0.0,
      "Point and unit tangent of a launching path at arc length s: straight, circular (plan radius > 0 turning "
      "left) or helical (a circle in plan on the grade dz/ds_plan), from `origin` with the plan heading [rad].");
  M.def(
      "launch_on_path",
      [launch_path](const PyModel& self, const std::vector<Id>& chain, const std::vector<double>& piers,
                    const std::vector<double>& advances, const std::string& kind, const Vec3T& origin, double heading,
                    double radius, double grade, int stations, double path_tol) {
        ModelPtr model = self.solver_model();
        LaunchOptions opt;
        opt.stations_per_bar = stations;
        opt.path_tol = path_tol;
        const LaunchPathResults r = engine([&] {
          return launch_on_path(*model, chain, launch_path(kind, origin, heading, radius, grade), piers, advances, opt);
        });
        py::list positions;
        for (const LaunchPosition& p : r.positions) {
          py::dict d;
          d["advance"] = p.advance;
          d["stable"] = p.stable;
          d["contact_s"] = p.contact_s;
          d["reaction"] = p.reaction;
          d["guide"] = p.guide;
          d["uplift"] = p.uplift;
          d["front_overhang"] = p.front_overhang;
          d["rear_overhang"] = p.rear_overhang;
          d["moment"] = p.moment;
          positions.append(d);
        }
        py::dict out;
        out["positions"] = positions;
        out["m_min"] = r.m_min;
        out["m_max"] = r.m_max;
        return out;
      },
      py::arg("chain"), py::arg("piers"), py::arg("advances"), py::arg("kind") = "straight",
      py::arg("origin") = Vec3T(0.0, 0.0, 0.0), py::arg("heading") = 0.0, py::arg("radius") = 0.0,
      py::arg("grade") = 0.0, py::arg("stations") = 5, py::arg("path_tol") = 1e-6,
      "Incremental launching of the girder `chain` (its bars from the rear end to the front end, nodes on the "
      "path) under its self-weight: piers at ground abscissae `piers`, one position per advance (the front end "
      "at ground abscissa `advance`). Each pier bears at its exact contact abscissa, between two nodes as well. "
      "Returns dict(positions=[dict(advance, stable, contact_s: contact abscissa along the girder from its "
      "rear end per pier (NaN off the girder), reaction: vertical reaction per pier, guide: lateral force per "
      "pier, moment: girder My at each contact, uplift: piers lifting off, front_overhang, rear_overhang)], "
      "m_min, m_max: envelopes of My over the stable positions, per bar station). The supports, cases and combinations of the model are ignored.");

  // ---------------- EN 1990 combinations ----------------
  M.def(
      "generate_combinations",
      [](PyModel& self, const py::list& permanent, const py::list& variable, std::string state, bool use_610ab,
         double xi) {
        std::vector<PermanentAction> G;
        for (const py::handle& h : permanent) {
          PermanentAction g;
          if (py::isinstance<py::dict>(h)) {
            const py::dict d = py::reinterpret_borrow<py::dict>(h);
            g.case_id = d["case"].cast<Id>();
            if (d.contains("gamma_sup")) g.gamma_sup = d["gamma_sup"].cast<double>();
            if (d.contains("gamma_inf")) g.gamma_inf = d["gamma_inf"].cast<double>();
          } else {
            g.case_id = h.cast<Id>();
          }
          G.push_back(g);
        }
        for (const PermanentAction& g : G)
          if (g.gamma_sup != G.front().gamma_sup || g.gamma_inf != G.front().gamma_inf)
            throw py::value_error("generate_combinations: the permanent actions take their partial factors together "
                                  "(all gamma_G,sup, then all gamma_G,inf): give them one pair");
        std::vector<VariableAction> Q;
        for (const py::handle& h : variable) Q.push_back(variable_of(h));
        for (char& c : state) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        std::vector<Id> ids = engine([&] {
          if (state == "uls") {
            Ec0Options o;
            o.use_610ab = use_610ab;
            o.xi = xi;
            return generate_uls(self.m, G, Q, o);
          }
          if (state == "characteristic") return generate_sls_characteristic(self.m, G, Q);
          if (state == "frequent") return generate_sls_frequent(self.m, G, Q);
          if (state == "quasi_permanent") return generate_sls_quasi_permanent(self.m, G, Q);
          throw Error("generate_combinations: state 'uls', 'characteristic', 'frequent' or 'quasi_permanent'");
        });
        for (const Combination& c : self.m.combos) self.next_combo = std::max(self.next_combo, c.id + 1);
        return ids;
      },
      py::arg("permanent"), py::arg("variable") = py::list(), py::arg("state") = "uls", py::arg("use_610ab") = false,
      py::arg("xi") = 0.85,
      "Adds the EN 1990 combinations to the model, each variable action leading in turn, and returns their ids. "
      "permanent = [case, ...] or [dict(case, gamma_sup=1.35, gamma_inf=1.0)] (one pair for all); variable = "
      "[dict(case, name, psi0, psi1, psi2, gamma_q)] or [(case, 'traffic' | 'wind' | 'thermal' | 'snow')] (French "
      "annex values). state 'uls' (6.10, or 6.10a + 6.10b with use_610ab and the factor xi), 'characteristic', "
      "'frequent', 'quasi_permanent'.");

  // ---------------- envelopes ----------------
  py::class_<PyBound> bound(mod, "Bound", "An envelope bound or a concomitant state: the accessors of a case result.");
  bind_case_accessors(bound, [](const PyBound& s) { return std::pair<ModelPtr, const CaseResult&>(s.model, *s.c); },
                      "end_forces");
  bound.def_property_readonly("label", [](const PyBound& s) { return s.label; });
  bound.def("__repr__", [](const PyBound& s) { return "<keta.Bound: " + s.label + ">"; });

  auto R = py::reinterpret_borrow<py::class_<PyResults>>(mod.attr("Results"));
  R.def(
      "envelope",
      [](const PyResults& s, std::optional<std::vector<Id>> cases, std::optional<std::vector<Id>> combos) {
        const auto src = sources_of(s, cases, combos, nullptr);
        auto env = std::make_shared<const Envelope>(engine([&] { return envelope(src); }));
        return std::make_pair(PyBound{s.model, std::shared_ptr<const CaseResult>(env, &env->lo.value), "minimum"},
                              PyBound{s.model, std::shared_ptr<const CaseResult>(env, &env->hi.value), "maximum"});
      },
      py::arg("cases") = py::none(), py::arg("combos") = py::none(),
      "Envelope of results (default: every case and combination): (minimum, maximum), each component taken "
      "independently (displacements, reactions, forces at every study section).");
  R.def(
      "envelope_factored",
      [](const PyResults& s, const std::vector<std::tuple<Id, double, double>>& terms) {
        std::vector<FactoredSource> src;
        for (const auto& [id, cmin, cmax] : terms) {
          const CaseResult* c = nullptr;
          for (const CaseResult& x : s.r->cases)
            if (x.id == id) c = &x;
          if (!c) throw py::key_error("no load case " + std::to_string(id));
          src.push_back(FactoredSource{c, cmin, cmax});
        }
        auto env = std::make_shared<const Envelope>(engine([&] { return envelope_factored(src); }));
        return std::make_pair(PyBound{s.model, std::shared_ptr<const CaseResult>(env, &env->lo.value), "minimum"},
                              PyBound{s.model, std::shared_ptr<const CaseResult>(env, &env->hi.value), "maximum"});
      },
      py::arg("terms"),
      "Envelope of factored cases: terms = [(case, c_min, c_max)]; for each component every case is taken ALONE "
      "with the factor of [c_min, c_max] that makes it largest (resp. smallest), and the bound is the extreme over "
      "the cases. It is not a combination (add_combination / generate_combinations sum the cases). Returns "
      "(minimum, maximum).");
  R.def(
      "concomitant",
      [](const PyResults& s, std::string what, Id bar, std::size_t station, int comp, Id node,
         std::optional<std::vector<Id>> cases, std::optional<std::vector<Id>> combos) {
        std::vector<std::string> names;
        const auto src = sources_of(s, cases, combos, &names);
        Privileged p;
        if (what == "bar_force") p.what = Privileged::What::StationForce;
        else if (what == "reaction") p.what = Privileged::What::Reaction;
        else if (what == "displacement") p.what = Privileged::What::Displacement;
        else throw py::value_error("concomitant: what 'bar_force', 'reaction' or 'displacement'");
        if (comp < 0 || comp > 5) throw py::value_error("concomitant: comp in 0..5");
        p.bar = bar;
        p.station = station;
        p.comp = comp;
        p.node = p.what == Privileged::What::Displacement ? static_cast<Id>(node_pos(*s.model, node)) : node;
        const Concomitant c = engine([&] { return envelope_concomitant(src, p); });
        auto keep = s.r;
        py::dict d;
        d["max"] = c.max;
        d["min"] = c.min;
        d["source_max"] = names[c.src_max];
        d["source_min"] = names[c.src_min];
        d["state_max"] = PyBound{s.model, std::shared_ptr<const CaseResult>(keep, c.state_max), names[c.src_max]};
        d["state_min"] = PyBound{s.model, std::shared_ptr<const CaseResult>(keep, c.state_min), names[c.src_min]};
        return d;
      },
      py::arg("what"), py::arg("bar") = 0, py::arg("station") = 0, py::arg("comp") = 0, py::arg("node") = 0,
      py::arg("cases") = py::none(), py::arg("combos") = py::none(),
      "Extremes of one component and the complete states that go with them: what = 'bar_force' (bar, station "
      "index, comp 0..5 = N Vy Vz T My Mz), 'reaction' or 'displacement' (node, comp 0..5). Returns dict(max, min, "
      "source_max, source_min, state_max, state_min) — the states give every other result of the source that "
      "governs.");

  // ---------------- EN 1998-1 spectra ----------------
  mod.def(
      "ec8_spectrum",
      [](double ag, double S, double TB, double TC, double TD, std::optional<double> q, double zeta, double beta,
         double plateau, std::optional<std::vector<double>> periods) {
        ResponseSpectrum sp = q ? ResponseSpectrum::ec8_design(ag, S, TB, TC, TD, *q)
                                : ResponseSpectrum::ec8_elastic(ag, S, TB, TC, TD);
        sp.beta_floor = beta;
        sp.plateau = plateau;
        std::vector<double> T;
        if (periods) {
          T = *periods;
        } else {
          for (int i = 0; i <= 400; ++i) T.push_back(0.01 * i);
          for (double c : {TB, TC, TD}) T.push_back(c);
          std::sort(T.begin(), T.end());
          T.erase(std::unique(T.begin(), T.end()), T.end());
        }
        std::vector<std::pair<double, double>> out;
        for (double t : T) out.emplace_back(t, engine([&] { return sp.sa(t, zeta); }));
        return out;
      },
      py::arg("ag"), py::arg("S") = 1.0, py::arg("TB") = 0.15, py::arg("TC") = 0.4, py::arg("TD") = 2.0,
      py::arg("q") = py::none(), py::arg("zeta") = 0.05, py::arg("beta") = 0.2, py::arg("plateau") = 2.5,
      py::arg("periods") = py::none(),
      "EN 1998-1 spectrum as a table [(T, Sa)] [s, m/s^2] (default periods 0 to 4 s by 0.01 s plus TB, TC, TD): "
      "elastic (eq. 3.2-3.5, damping correction eta(zeta)) or, with the behaviour factor q, design (eq. 3.13-3.16, "
      "floor beta ag). plateau = 3.0 for a vertical spectrum. The table is the `spectrum` argument of "
      "response_spectrum.");

  // ---------------- EN 1993-1-1 / 1-5 values ----------------
  auto curve_of = [](const std::string& c) {
    if (c == "a0") return ec3::Curve::a0;
    if (c == "a") return ec3::Curve::a;
    if (c == "b") return ec3::Curve::b;
    if (c == "c") return ec3::Curve::c;
    if (c == "d") return ec3::Curve::d;
    throw py::value_error("buckling curve 'a0', 'a', 'b', 'c' or 'd'");
  };
  mod.def("ec3_epsilon", [](double fy) { return engine([&] { return ec3::epsilon(fy); }); }, py::arg("fy"),
          "EN 1993-1-1 Table 5.2: epsilon = sqrt(235 MPa / fy), fy [Pa].");
  mod.def(
      "ec3_class",
      [](double c_t, double fy, std::string part, double alpha, double psi, bool tip_in_compression) {
        const double eps = ec3::epsilon(fy);
        if (part == "internal") return engine([&] { return ec3::class_internal(c_t, eps, alpha, psi); });
        if (part == "outstand")
          return engine([&] { return ec3::class_outstand(c_t, eps, alpha, psi, tip_in_compression); });
        throw py::value_error("ec3_class: part 'internal' or 'outstand'");
      },
      py::arg("c_t"), py::arg("fy"), py::arg("part") = "internal", py::arg("alpha") = 1.0, py::arg("psi") = 1.0,
      py::arg("tip_in_compression") = true,
      "Class (1 to 4) of a compression part, EN 1993-1-1 Table 5.2: slenderness c/t, fy [Pa]; alpha = compressed "
      "fraction under the plastic distribution, psi = sigma2/sigma1 of the elastic one (1, 1: pure compression).");
  mod.def(
      "ec3_plate_reduction",
      [](double b_t, double fy, double psi, std::string part, bool max_compression_at_free_edge) {
        const double eps = ec3::epsilon(fy);
        py::dict d;
        double ks = 0, lp = 0, rho = 0;
        engine([&] {
          if (part == "internal") {
            ks = ec3::ksigma_internal(psi);
            lp = ec3::lambda_p(b_t, eps, ks);
            rho = ec3::rho_internal(lp, psi);
          } else if (part == "outstand") {
            ks = ec3::ksigma_outstand(psi, max_compression_at_free_edge);
            lp = ec3::lambda_p(b_t, eps, ks);
            rho = ec3::rho_outstand(lp);
          } else {
            throw Error("ec3_plate_reduction: part 'internal' or 'outstand'");
          }
          return 0;
        });
        d["k_sigma"] = ks;
        d["lambda_p"] = lp;
        d["rho"] = rho;
        return d;
      },
      py::arg("b_t"), py::arg("fy"), py::arg("psi") = 1.0, py::arg("part") = "internal",
      py::arg("max_compression_at_free_edge") = true,
      "EN 1993-1-5 §4.4 effective width of a plate: dict(k_sigma (Tables 4.1 / 4.2), lambda_p, rho (eq. 4.2 / "
      "4.3)).");
  mod.def(
      "ec3_chi",
      [curve_of](double lambda_bar, const std::string& curve, bool lateral_torsional_rolled, double lambda0,
                 double beta) {
        const double a = ec3::imperfection(curve_of(curve));
        return engine([&] {
          return lateral_torsional_rolled ? ec3::chi_lt_rolled(lambda_bar, a, lambda0, beta) : ec3::chi(lambda_bar, a);
        });
      },
      py::arg("lambda_bar"), py::arg("curve"), py::arg("lateral_torsional_rolled") = false, py::arg("lambda0") = 0.4,
      py::arg("beta") = 0.75,
      "Buckling reduction factor chi of EN 1993-1-1: eq. 6.49 (flexural, and lateral-torsional general case) for "
      "the curve 'a0' to 'd' (imperfection 0.13 to 0.76, Tables 6.1 / 6.3); lateral_torsional_rolled: eq. 6.57 "
      "with lambda0 and beta.");
  mod.def(
      "ec3_buckling_curve",
      [](double h, double b, double tf, bool rolled, std::string kind, bool about_y) {
        ec3::Curve c;
        if (kind == "flexural") c = ec3::curve_flexural_i(rolled, h, b, tf, about_y);
        else if (kind == "lt_general") c = ec3::curve_lt_general_i(rolled, h, b);
        else if (kind == "lt_rolled") c = ec3::curve_lt_rolled_i(rolled, h, b);
        else throw py::value_error("ec3_buckling_curve: kind 'flexural', 'lt_general' or 'lt_rolled'");
        static const char* names[] = {"a0", "a", "b", "c", "d"};
        return std::string(names[static_cast<int>(c)]);
      },
      py::arg("h"), py::arg("b"), py::arg("tf"), py::arg("rolled") = true, py::arg("kind") = "flexural",
      py::arg("about_y") = true,
      "Buckling curve of an I section, EN 1993-1-1 Table 6.2 (flexural), 6.4 (lateral-torsional, general) or 6.5 "
      "(lateral-torsional, rolled or equivalent welded); dimensions [m].");
}
