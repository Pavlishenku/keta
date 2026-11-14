# Beam fixed at A, simply supported at B, cantilevering from B to C, force P at the tip C.
import keta

E, A, I = 205e9, 81.9e-4, 22964.9e-8     # 400 mm H-section, steel
L, P = 4.0, 10e3                          # length AB = BC, downward force at C

m = keta.Model()
steel = m.add_material(e=E, nu=0.3)
sec = m.add_section_props(sx=A, iy=I, iz=I, ix=I)
a, b, c = m.add_node(0, 0, 0), m.add_node(L, 0, 0), m.add_node(2 * L, 0, 0)
ab, bc = m.add_beam(a, b, sec, steel), m.add_beam(b, c, sec, steel)   # one bar per span
m.fix(a)
m.add_support(b, dz=True)
for n in (b, c):
    m.add_support(n, dy=True, rx=True, rz=True)    # plane problem in XZ
m.add_point_load(c, fz=-P)
case = m.analyze().case(1)

R = dict(case.reactions())
values = [  # (quantity, keta, closed form)
    ("My at A [kN.m]", case.bar_forces(ab)[0][5] / 1e3, -P * L / 2 / 1e3),
    ("My at B [kN.m]", case.bar_forces(bc)[0][5] / 1e3, P * L / 1e3),
    ("Vz in AB [kN]", case.bar_forces(ab)[0][3] / 1e3, 1.5 * P / 1e3),
    ("reaction A [kN]", R[a][2] / 1e3, -1.5 * P / 1e3),
    ("reaction B [kN]", R[b][2] / 1e3, 2.5 * P / 1e3),
    ("deflection C [mm]", case.displacement(c)[2] * 1e3, -7 * P * L**3 / (12 * E * I) * 1e3),
]
for name, v, ref in values:
    diff = 100 * (v - ref) / ref
    print(f"{name + ':':18s} keta {v:9.4f}   closed form {ref:9.4f}   diff. {round(diff, 4) + 0:.4f} %")
    assert abs(diff) < 1e-7
