# Charges sur barres : concentrée, répartie partielle (trapèze) et cas de charge nommés.
import keta

L = 8.0
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
poutre = m.add_beam(a, b, m.add_section("HEB300"), m.steel())
m.add_support(a, dx=True, dy=True, dz=True, rx=True)
m.add_support(b, dy=True, dz=True)

m.add_bar_point_load(poutre, at=0.25, fz=-40e3, case=1)                 # 40 kN au quart
m.add_linear_load(poutre, w1=(0, 0, 0), w2=(0, 0, -30e3), x1=0.5, x2=1.0, case=2)   # triangle
m.add_udl(poutre, wz=-5e3, case=3)
m.add_point_load(b, my=20e3, case=4)                                    # couple sur appui
for c, nom in [(1, "concentrée"), (2, "triangle"), (3, "uniforme"), (4, "couple")]:
    m.name_case(c, nom)

r = m.analyze(stations=17)
print("cas :", r.case_ids, "-", r.case(1))
for c in r.case_ids:
    cas = r.case(c)
    Ra, Rb = (rz[2] for _, rz in cas.reactions())
    Mmax = max(cas.bar_forces(poutre), key=lambda s: abs(s[5]))
    print(f"{c} RA = {Ra / 1e3:7.3f} kN  RB = {Rb / 1e3:7.3f} kN  "
          f"|M|max = {abs(Mmax[5]) / 1e3:7.3f} kN.m à x/L = {Mmax[0]:.4f}")
# Contrôles : RA(1) = 3P/4 ; triangle de 0,5L à L : R = 30 x 4 / 2 = 60 kN à 2/3 de [4, 8] m
assert abs(dict(r.case(1).reactions())[a][2] - 30e3) < 1e-6
assert abs(dict(r.case(2).reactions())[b][2] - 60e3 * (4 + 8 / 3) / 8) < 1e-6
