# Poutre encastrée en A, appuyée en B, en console de B à C, force P à l'extrémité C.
import keta

E, A, I = 205e9, 81.9e-4, 22964.9e-8     # profilé en H de 400 mm, acier
L, P = 4.0, 10e3                          # longueur AB = BC, force descendante en C

m = keta.Model()
acier = m.add_material(e=E, nu=0.3)
sec = m.add_section_props(sx=A, iy=I, iz=I, ix=I)
a, b, c = m.add_node(0, 0, 0), m.add_node(L, 0, 0), m.add_node(2 * L, 0, 0)
ab, bc = m.add_beam(a, b, sec, acier), m.add_beam(b, c, sec, acier)   # une barre par travée
m.fix(a)
m.add_support(b, dz=True)
for n in (b, c):
    m.add_support(n, dy=True, rx=True, rz=True)    # problème plan XZ
m.add_point_load(c, fz=-P)
cas = m.analyze().case(1)

R = dict(cas.reactions())
valeurs = [  # (grandeur, keta, forme fermée)
    ("My en A  [kN.m]", cas.bar_forces(ab)[0][5] / 1e3, -P * L / 2 / 1e3),
    ("My en B  [kN.m]", cas.bar_forces(bc)[0][5] / 1e3, P * L / 1e3),
    ("Vz sur AB  [kN]", cas.bar_forces(ab)[0][3] / 1e3, 1.5 * P / 1e3),
    ("réaction A [kN]", R[a][2] / 1e3, -1.5 * P / 1e3),
    ("réaction B [kN]", R[b][2] / 1e3, 2.5 * P / 1e3),
    ("flèche C   [mm]", cas.displacement(c)[2] * 1e3, -7 * P * L**3 / (12 * E * I) * 1e3),
]
for nom, v, ref in valeurs:
    ecart = 100 * (v - ref) / ref
    print(f"{nom} : keta {v:9.4f}   forme fermée {ref:9.4f}   écart {round(ecart, 4) + 0:.4f} %")
    assert abs(ecart) < 1e-7
