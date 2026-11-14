# Portique deux fois hyperstatique : montant articulé en pied (1), montant court encastré (2).
# Comparaison avec la méthode des forces (inconnues : réactions X1 horizontale, X2 verticale en 1).
import keta

E, I, q, F = 30e9, 1e-3, 1e3, 3e3     # q horizontale sur le montant 1-3, F au milieu de la traverse
m = keta.Model()
beton = m.add_material(e=E, nu=0.2)
double = m.add_section_props(sx=1e3, iy=2 * I, iz=2 * I, ix=I)   # aire très grande :
simple = m.add_section_props(sx=1e3, iy=I, iz=I, ix=I)           # barres axialement rigides
n1, n3, c, n4, n2 = (m.add_node(x, 0, z) for x, z in ((0, 0), (0, 6), (3, 6), (6, 6), (6, 3)))
p13 = m.add_beam(n1, n3, double, beton)
t3c, tc4 = m.add_beam(n3, c, double, beton), m.add_beam(c, n4, double, beton)
p42 = m.add_beam(n4, n2, simple, beton)
m.pin(n1)
m.fix(n2)
for n in (n1, n3, c, n4):
    m.add_support(n, dy=True, rx=True, rz=True)   # portique plan XZ
m.add_udl(p13, wx=q)
m.add_point_load(c, fz=-F)
cas = m.analyze(stations=601).case(1)            # 601 sections d'étude : pas de 1 cm

# Méthode des forces, coefficients multipliés par EI (kN, m) : d11 X1 + d12 X2 + D1 = 0, etc.
d11, d22, d12, D1, D2 = 207.0, 144.0, -135.0, 702.0, -519.75
det = d11 * d22 - d12**2
X1, X2 = (-D1 * d22 + d12 * D2) / det, (-d11 * D2 + d12 * D1) / det       # kN
# M = X1 M1 + X2 M2 + M0 aux points 3, C, 4 et 2 ; maximum du montant 1-3 en z = -X1
ref = {"X1 [kN]": X1, "X2 [kN]": X2, "M3 [kN.m]": 6 * X1 + 18, "MC [kN.m]": 6 * X1 - 3 * X2 + 18,
       "M4 [kN.m]": 6 * X1 - 6 * X2 + 27, "M2 [kN.m]": 3 * X1 - 6 * X2 + 9,
       "Mmax [kN.m]": X1**2 / 2}

R1 = dict(cas.reactions())[n1]
Mp = [abs(s[5]) for s in cas.bar_forces(p13)]
val = {"X1 [kN]": R1[0], "X2 [kN]": R1[2], "M3 [kN.m]": cas.bar_forces(p13)[-1][5],
       "MC [kN.m]": cas.bar_forces(t3c)[-1][5], "M4 [kN.m]": cas.bar_forces(tc4)[-1][5],
       "M2 [kN.m]": cas.bar_forces(p42)[-1][5], "Mmax [kN.m]": max(Mp)}
for nom, ref_v in ref.items():
    v = val[nom] / 1e3
    if nom.startswith("M"):                       # moments comparés en valeur absolue
        v, ref_v = abs(v), abs(ref_v)
    ecart = 100 * (v - ref_v) / ref_v
    print(f"{nom:12s} keta {v:8.4f}  méthode des forces {ref_v:8.4f}  "
          f"écart {round(ecart, 4) + 0:.4f} %")
    assert abs(ecart) < 1e-3
z = Mp.index(max(Mp)) * 0.01
print(f"moment maximal du montant 1-3 en z = {z:.2f} m (forme fermée : z = -X1 = {-X1:.4f} m)")
