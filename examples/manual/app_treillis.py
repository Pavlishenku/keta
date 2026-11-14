# Statically determinate plane truss of four pin-jointed bars: axial forces and displacements.
import math
import keta

E, F = 1.962e11, 9.81e3                    # Young's modulus, downward vertical force at D
P = {"A": (0.0, 0.0), "B": (1.0, 0.0), "C": (0.5, 0.5), "D": (2.0, 1.0)}
BARS = {"AC": 2e-4, "BC": 2e-4, "CD": 1e-4, "BD": 1e-4}         # areas [m2]

m = keta.Model()
steel = m.add_material(e=E, nu=0.3)
nd = {k: m.add_node(x, 0, z) for k, (x, z) in P.items()}
bar = {}
for name, area in BARS.items():
    sec = m.add_section_props(sx=area, iy=1e-12, iz=1e-12, ix=1e-12)
    bar[name] = m.add_truss(nd[name[0]], nd[name[1]], sec, steel)
for k in "AB":
    m.add_support(nd[k], dx=True, dy=True, dz=True)
for k in "CD":
    m.add_support(nd[k], dy=True)                                 # plane truss in XZ
m.add_point_load(nd["D"], fz=-F)
case = m.analyze().case(1)

# Closed form: equilibrium of nodes D then C; compatibility of the elongations e = N L / (E A)
r2, r10 = math.sqrt(2), math.sqrt(10)
N = {"AC": r2 * F, "BC": -F / r2, "CD": r10 / 2 * F, "BD": -3 / r2 * F}
def unit_vector(name):
    (x1, z1), (x2, z2) = P[name[0]], P[name[1]]
    L = math.hypot(x2 - x1, z2 - z1)
    return (x2 - x1) / L, (z2 - z1) / L, L
def displacement(b1, e1, b2, e2):                                 # n1 . u = e1 and n2 . u = e2
    (a1, c1, _), (a2, c2, _) = unit_vector(b1), unit_vector(b2)
    det = a1 * c2 - a2 * c1
    return (e1 * c2 - e2 * c1) / det, (a1 * e2 - a2 * e1) / det
e = {k: N[k] * unit_vector(k)[2] / (E * BARS[k]) for k in BARS}
uC = displacement("AC", e["AC"], "BC", e["BC"])
nx, nz, _ = unit_vector("CD")
uD = displacement("BD", e["BD"], "CD", e["CD"] + nx * uC[0] + nz * uC[1])

for k in BARS:
    v = case.bar_forces(bar[k])[0][1]
    print(f"N {k}: keta {v / 1e3:9.4f} kN   closed form {N[k] / 1e3:9.4f} kN   "
          f"diff. {round(100 * (v - N[k]) / N[k], 4) + 0:.4f} %")
    assert abs(v - N[k]) < 1e-9 * abs(N[k])
for k, u in (("C", uC), ("D", uD)):
    ux, _, uz, *_ = case.displacement(nd[k])
    for name, v, ref in ((f"u{k}", ux, u[0]), (f"w{k}", uz, u[1])):
        print(f"{name}: keta {v * 1e3:9.6f} mm   closed form {ref * 1e3:9.6f} mm   "
              f"diff. {round(100 * (v - ref) / ref, 4) + 0:.4f} %")
        assert abs(v - ref) < 1e-9 * abs(ref)
