# Suspension bridge: a deck of 21 spans of 10 m hung from a main cable between two saddles
# 240 m apart. Form finding of the cable under the hanger loads, assembly, then
# hangers at their unstressed length, two-stage correction of the unstressed lengths of the cable (the stiff
# deck carries part of the load).
import keta

def bridge():
    m = keta.Model()
    steel = m.add_material(200e9, 0.3, rho=0.0)
    wire = m.add_material(195e9, 0.3, rho=7850.0)
    rod = m.add_material(160e9, 0.3, rho=0.0)
    s_deck = m.add_section_props(0.5, 0.8, 5.0, 1.0)
    s_cable = m.add_section_props(0.30, 1e-6, 1e-6, 1e-6)
    s_hanger = m.add_section_props(0.01, 1e-8, 1e-8, 1e-8)
    nh, x0, dx, q = 21, 20.0, 10.0, 150e3
    nodes = [m.add_node(x0 - dx / 2, 0, 5.0)] + [m.add_node(x0 + dx * i, 0, 5.0) for i in range(nh)] \
        + [m.add_node(x0 + dx * (nh - 1) + dx / 2, 0, 5.0)]
    bars = [m.add_beam(a, b, s_deck, steel) for a, b in zip(nodes, nodes[1:])]
    m.add_support(nodes[0], dx=True, dy=True, rx=True)       # the deck is hung, not supported
    m.add_support(nodes[-1], dy=True)
    for b in bars:
        m.add_udl(b, wz=-q, case=1)
    data = dict(left=(0, 0, 40), right=(240, 0, 40), hanger_x=[x0 + dx * i for i in range(nh)],
                hanger_load=[q * dx] * nh, hanger_deck_nodes=nodes[1:-1], cable_material=wire,
                cable_section=s_cable, hanger_material=rod, hanger_section=s_hanger, target_sag=25.0)
    return m, data, q * dx * nh

# 1. form finding alone: horizontal tension and unstressed length of the cable
m, data, weight = bridge()
shape = keta.form_find_cable(data["left"], data["right"], data["hanger_x"], data["hanger_load"],
                             ea=195e9 * 0.30, w=7850.0 * 0.30 * 9.80665, target_sag=25.0)
print(f"shape: H = {shape['H'] / 1e6:.2f} MN, sag {shape['sag']:.3f} m, unstressed length {shape['total_l0']:.3f} m")

# 2. assembly in the model and equilibrium under the deck weight
sm = m.build_suspension(**data)
r = m.cable_static(case=1)
deviation = max(abs(r.displacement(n)[2]) for n in sm["cable_nodes"][1:-1])
print(f"assembled: {len(sm['cables'])} cable segments, {len(sm['hanger_bars'])} hangers; "
      f"vertical reaction {r.sum_reactions()[2] / 1e6:.3f} MN for {weight / 1e6:.3f} MN of deck; "
      f"cable off its profile by {deviation * 1e3:.2f} mm")

# 3. the same hangers at their unstressed length: their length defect goes into the permanent case
m1, data1, _ = bridge()
sm1 = m1.build_suspension(**data1, hanger_fit_case=1)
r1 = m1.cable_static(case=1)
defl = min(r.displacement(n)[2] for n in data["hanger_deck_nodes"])
defl1 = min(r1.displacement(n)[2] for n in data1["hanger_deck_nodes"])
print(f"hangers: placed stress-free, the deck goes down by {-defl * 1e3:.1f} mm; at their unstressed length "
      f"(from {min(sm1['hanger_l0']):.3f} to {max(sm1['hanger_l0']):.3f} m), by {-defl1 * 1e3:.1f} mm")

# 4. two-stage correction: the unstressed lengths bring the cable back onto its profile
m2, data2, _ = bridge()
u = m2.solve_uelm(1, **data2, tol=1e-10)
r2 = m2.cable_static(case=1)
deviation2 = max(abs(r2.displacement(n)[2]) for n in u["cable_nodes"][1:-1])
print(f"two stages: deviation {u['deviation_stage1'] * 1e3:.2f} mm before correction, {deviation2 * 1e3:.4f} mm after, "
      f"{u['iterations']} iterations; total unstressed length {sum(u['l0']):.3f} m")
assert deviation2 < 1e-7 and abs(r.sum_reactions()[2] - weight - 7850.0 * 0.30 * 9.80665 * shape["total_l0"]) < 1e-6 * weight
