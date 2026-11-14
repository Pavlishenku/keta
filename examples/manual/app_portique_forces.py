# Portal frame statically indeterminate to the second degree: post pinned at its foot (1), short
# post fixed (2). Comparison with the force method (unknowns: reactions X1, horizontal, and X2,
# vertical, at 1).
import keta

E, I, q, F = 30e9, 1e-3, 1e3, 3e3     # q horizontal on post 1-3, F at mid-span of the girder
m = keta.Model()
concrete = m.add_material(e=E, nu=0.2)
double = m.add_section_props(sx=1e3, iy=2 * I, iz=2 * I, ix=I)   # very large area:
single = m.add_section_props(sx=1e3, iy=I, iz=I, ix=I)           # axially rigid bars
n1, n3, c, n4, n2 = (m.add_node(x, 0, z) for x, z in ((0, 0), (0, 6), (3, 6), (6, 6), (6, 3)))
p13 = m.add_beam(n1, n3, double, concrete)
g3c, gc4 = m.add_beam(n3, c, double, concrete), m.add_beam(c, n4, double, concrete)
p42 = m.add_beam(n4, n2, single, concrete)
m.pin(n1)
m.fix(n2)
for n in (n1, n3, c, n4):
    m.add_support(n, dy=True, rx=True, rz=True)   # plane frame in XZ
m.add_udl(p13, wx=q)
m.add_point_load(c, fz=-F)
case = m.analyze(stations=601).case(1)           # 601 stations: 1 cm spacing

# Force method, coefficients multiplied by EI (kN, m): d11 X1 + d12 X2 + D1 = 0, etc.
d11, d22, d12, D1, D2 = 207.0, 144.0, -135.0, 702.0, -519.75
det = d11 * d22 - d12**2
X1, X2 = (-D1 * d22 + d12 * D2) / det, (-d11 * D2 + d12 * D1) / det       # kN
# M = X1 M1 + X2 M2 + M0 at points 3, C, 4 and 2; maximum of post 1-3 at z = -X1
ref = {"X1 [kN]": X1, "X2 [kN]": X2, "M3 [kN.m]": 6 * X1 + 18, "MC [kN.m]": 6 * X1 - 3 * X2 + 18,
       "M4 [kN.m]": 6 * X1 - 6 * X2 + 27, "M2 [kN.m]": 3 * X1 - 6 * X2 + 9,
       "Mmax [kN.m]": X1**2 / 2}

R1 = dict(case.reactions())[n1]
Mp = [abs(s[5]) for s in case.bar_forces(p13)]
val = {"X1 [kN]": R1[0], "X2 [kN]": R1[2], "M3 [kN.m]": case.bar_forces(p13)[-1][5],
       "MC [kN.m]": case.bar_forces(g3c)[-1][5], "M4 [kN.m]": case.bar_forces(gc4)[-1][5],
       "M2 [kN.m]": case.bar_forces(p42)[-1][5], "Mmax [kN.m]": max(Mp)}
for name, ref_v in ref.items():
    v = val[name] / 1e3
    if name.startswith("M"):                      # moments compared in absolute value
        v, ref_v = abs(v), abs(ref_v)
    diff = 100 * (v - ref_v) / ref_v
    print(f"{name:12s} keta {v:8.4f}  force method {ref_v:8.4f}  diff. {round(diff, 4) + 0:.4f} %")
    assert abs(diff) < 1e-3
z = Mp.index(max(Mp)) * 0.01
print(f"maximum moment of post 1-3 at z = {z:.2f} m (closed form: z = -X1 = {-X1:.4f} m)")
