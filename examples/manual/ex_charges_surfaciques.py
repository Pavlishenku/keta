# Area loads: load distribution zone on a panel of beams, polygonal
# pressure and hydrostatic pressure on shells.
import keta

# 1) Floor of 6 m x 4 m bordered by four beams: 5 kN/m2 distributed to the beams
for mode in ("two_way", "one_way"):
    m = keta.Model()
    c = [m.add_node(x, y, 0) for x, y in ((0, 0), (6, 0), (6, 4), (0, 4))]
    sec, mat = m.add_section("IPE300"), m.steel()
    beams = [m.add_beam(c[i], c[(i + 1) % 4], sec, mat) for i in range(4)]
    for nd in c:
        m.pin(nd)
    nb = m.add_area_load([(0, 0, 0), (6, 0, 0), (6, 4, 0), (0, 4, 0)], qz=-5e3, mode=mode, span=(0, 1, 0))
    case_ = m.analyze().case(1)
    M = [abs(case_.bar_forces(p)[5][5]) / 1e3 for p in beams]
    print(f"{mode:8}: {nb} bar loads, total = {case_.sum_reactions()[2] / 1e3:.2f} kN, "
          f"midspan M of the beams = {', '.join(f'{x:.2f}' for x in M)} kN.m")

# 2) Vertical wall of 3 m x 2.5 m (shells), clamped at the base, retaining 2.5 m of water
m = keta.Model()
concrete = m.add_material(e=30e9, nu=0.2)
g = [[m.add_node(0.5 * i, 0.0, 0.5 * j) for i in range(7)] for j in range(6)]
for j in range(5):
    for i in range(6):
        m.add_shell(g[j][i], g[j][i + 1], g[j + 1][i + 1], 0.25, concrete, slab=1)
        m.add_shell(g[j][i], g[j + 1][i + 1], g[j + 1][i], 0.25, concrete, slab=1)
for nd in g[0]:
    m.fix(nd)
m.add_hydrostatic(gamma=9810.0, surface_point=(0.0, 0.0, 2.5), slab=1, case=1)
Fy = m.analyze().case(1).sum_reactions()[1]
print(f"hydrostatic thrust: |Fy| = {abs(Fy) / 1e3:.3f} kN (gamma h² b / 2 = {9810 * 2.5**2 * 3 / 2 / 1e3:.3f} kN)")

# 3) Pressure on a polygon with an opening, applied in plan on a horizontal slab
m = keta.Model()
concrete = m.add_material(e=30e9, nu=0.2)
g = [[m.add_node(1.0 * i, 1.0 * j, 0.0) for i in range(7)] for j in range(5)]
for j in range(4):
    for i in range(6):
        m.add_shell(g[j][i], g[j][i + 1], g[j + 1][i + 1], 0.2, concrete, slab=1)
        m.add_shell(g[j][i], g[j + 1][i + 1], g[j + 1][i], 0.2, concrete, slab=1)
for nd in g[0] + g[-1]:
    m.pin(nd)
outline = [(0.5, 0.5, 0), (5.5, 0.5, 0), (5.5, 3.5, 0), (0.5, 3.5, 0)]
opening = [(2.0, 1.5, 0), (3.0, 1.5, 0), (3.0, 2.5, 0), (2.0, 2.5, 0)]
m.add_polygon_pressure(outline, qz=-4e3, holes=[opening], case=1)
Fz = m.analyze().case(1).sum_reactions()[2]
print(f"polygonal pressure: Fz = {Fz / 1e3:.3f} kN (4 kN/m2 x (15 - 1) m2 = 56.000 kN)")
