# Treillis plan isostatique de quatre barres articulées : efforts normaux et déplacements.
import math
import keta

E, F = 1.962e11, 9.81e3                    # module d'Young, force verticale descendante en D
P = {"A": (0.0, 0.0), "B": (1.0, 0.0), "C": (0.5, 0.5), "D": (2.0, 1.0)}
BARRES = {"AC": 2e-4, "BC": 2e-4, "CD": 1e-4, "BD": 1e-4}       # aires [m2]

m = keta.Model()
acier = m.add_material(e=E, nu=0.3)
nd = {k: m.add_node(x, 0, z) for k, (x, z) in P.items()}
barre = {}
for nom, aire in BARRES.items():
    sec = m.add_section_props(sx=aire, iy=1e-12, iz=1e-12, ix=1e-12)
    barre[nom] = m.add_truss(nd[nom[0]], nd[nom[1]], sec, acier)
for k in "AB":
    m.add_support(nd[k], dx=True, dy=True, dz=True)
for k in "CD":
    m.add_support(nd[k], dy=True)                                 # treillis plan XZ
m.add_point_load(nd["D"], fz=-F)
cas = m.analyze().case(1)

# Forme fermée : équilibre des nœuds D puis C ; compatibilité des allongements e = N L / (E A)
r2, r10 = math.sqrt(2), math.sqrt(10)
N = {"AC": r2 * F, "BC": -F / r2, "CD": r10 / 2 * F, "BD": -3 / r2 * F}
def unitaire(nom):
    (x1, z1), (x2, z2) = P[nom[0]], P[nom[1]]
    L = math.hypot(x2 - x1, z2 - z1)
    return (x2 - x1) / L, (z2 - z1) / L, L
def deplacement(b1, e1, b2, e2):                                  # n1 . u = e1 et n2 . u = e2
    (a1, c1, _), (a2, c2, _) = unitaire(b1), unitaire(b2)
    det = a1 * c2 - a2 * c1
    return (e1 * c2 - e2 * c1) / det, (a1 * e2 - a2 * e1) / det
e = {k: N[k] * unitaire(k)[2] / (E * BARRES[k]) for k in BARRES}
uC = deplacement("AC", e["AC"], "BC", e["BC"])
nx, nz, _ = unitaire("CD")
uD = deplacement("BD", e["BD"], "CD", e["CD"] + nx * uC[0] + nz * uC[1])

for k in BARRES:
    v = cas.bar_forces(barre[k])[0][1]
    print(f"N {k} : keta {v / 1e3:9.4f} kN   forme fermée {N[k] / 1e3:9.4f} kN   "
          f"écart {round(100 * (v - N[k]) / N[k], 4) + 0:.4f} %")
    assert abs(v - N[k]) < 1e-9 * abs(N[k])
for k, u in (("C", uC), ("D", uD)):
    ux, _, uz, *_ = cas.displacement(nd[k])
    for nom, v, ref in ((f"u{k}", ux, u[0]), (f"w{k}", uz, u[1])):
        print(f"{nom} : keta {v * 1e3:9.6f} mm   forme fermée {ref * 1e3:9.6f} mm   "
              f"écart {round(100 * (v - ref) / ref, 4) + 0:.4f} %")
        assert abs(v - ref) < 1e-9 * abs(ref)
