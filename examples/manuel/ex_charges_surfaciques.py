# Charges surfaciques : zone de répartition sur un panneau de poutres, pression
# polygonale et pression hydrostatique sur des coques.
import keta

# 1) Plancher de 6 m x 4 m bordé de quatre poutres : 5 kN/m2 répartis sur les poutres
for mode in ("two_way", "one_way"):
    m = keta.Model()
    c = [m.add_node(x, y, 0) for x, y in ((0, 0), (6, 0), (6, 4), (0, 4))]
    sec, mat = m.add_section("IPE300"), m.steel()
    poutres = [m.add_beam(c[i], c[(i + 1) % 4], sec, mat) for i in range(4)]
    for nd in c:
        m.pin(nd)
    nb = m.add_area_load([(0, 0, 0), (6, 0, 0), (6, 4, 0), (0, 4, 0)], qz=-5e3, mode=mode, span=(0, 1, 0))
    cas = m.analyze().case(1)
    M = [abs(cas.bar_forces(p)[5][5]) / 1e3 for p in poutres]
    print(f"{mode:8}: {nb} charges de barre, total = {cas.sum_reactions()[2] / 1e3:.2f} kN, "
          f"M à mi-portée des poutres = {', '.join(f'{x:.2f}' for x in M)} kN.m")

# 2) Voile vertical de 3 m x 2,5 m (coques), encastré en pied, retenant 2,5 m d'eau
m = keta.Model()
beton = m.add_material(e=30e9, nu=0.2)
g = [[m.add_node(0.5 * i, 0.0, 0.5 * j) for i in range(7)] for j in range(6)]
for j in range(5):
    for i in range(6):
        m.add_shell(g[j][i], g[j][i + 1], g[j + 1][i + 1], 0.25, beton, slab=1)
        m.add_shell(g[j][i], g[j + 1][i + 1], g[j + 1][i], 0.25, beton, slab=1)
for nd in g[0]:
    m.fix(nd)
m.add_hydrostatic(gamma=9810.0, surface_point=(0.0, 0.0, 2.5), slab=1, case=1)
Fy = m.analyze().case(1).sum_reactions()[1]
print(f"poussée hydrostatique : |Fy| = {abs(Fy) / 1e3:.3f} kN (gamma h² b / 2 = {9810 * 2.5**2 * 3 / 2 / 1e3:.3f} kN)")

# 3) Pression sur un polygone avec trémie, appliquée en plan sur une dalle horizontale
m = keta.Model()
beton = m.add_material(e=30e9, nu=0.2)
g = [[m.add_node(1.0 * i, 1.0 * j, 0.0) for i in range(7)] for j in range(5)]
for j in range(4):
    for i in range(6):
        m.add_shell(g[j][i], g[j][i + 1], g[j + 1][i + 1], 0.2, beton, slab=1)
        m.add_shell(g[j][i], g[j + 1][i + 1], g[j + 1][i], 0.2, beton, slab=1)
for nd in g[0] + g[-1]:
    m.pin(nd)
contour = [(0.5, 0.5, 0), (5.5, 0.5, 0), (5.5, 3.5, 0), (0.5, 3.5, 0)]
tremie = [(2.0, 1.5, 0), (3.0, 1.5, 0), (3.0, 2.5, 0), (2.0, 2.5, 0)]
m.add_polygon_pressure(contour, qz=-4e3, holes=[tremie], case=1)
Fz = m.analyze().case(1).sum_reactions()[2]
print(f"pression polygonale : Fz = {Fz / 1e3:.3f} kN (4 kN/m2 x (15 - 1) m2 = 56.000 kN)")
