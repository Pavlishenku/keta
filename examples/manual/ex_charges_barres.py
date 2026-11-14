# Loads on bars: point, partial distributed (trapezoid) and named load cases.
import keta

L = 8.0
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
beam = m.add_beam(a, b, m.add_section("HEB300"), m.steel())
m.add_support(a, dx=True, dy=True, dz=True, rx=True)
m.add_support(b, dy=True, dz=True)

m.add_bar_point_load(beam, at=0.25, fz=-40e3, case=1)                   # 40 kN at quarter span
m.add_linear_load(beam, w1=(0, 0, 0), w2=(0, 0, -30e3), x1=0.5, x2=1.0, case=2)   # triangle
m.add_udl(beam, wz=-5e3, case=3)
m.add_point_load(b, my=20e3, case=4)                                    # moment on a support
for c, name in [(1, "point"), (2, "triangle"), (3, "uniform"), (4, "moment")]:
    m.name_case(c, name)

r = m.analyze(stations=17)
print("cases:", r.case_ids, "-", r.case(1))
for c in r.case_ids:
    case_ = r.case(c)
    Ra, Rb = (rz[2] for _, rz in case_.reactions())
    Mmax = max(case_.bar_forces(beam), key=lambda s: abs(s[5]))
    print(f"{c} RA = {Ra / 1e3:7.3f} kN  RB = {Rb / 1e3:7.3f} kN  "
          f"|M|max = {abs(Mmax[5]) / 1e3:7.3f} kN.m at x/L = {Mmax[0]:.4f}")
# Checks: RA(1) = 3P/4; triangle from 0.5L to L: R = 30 x 4 / 2 = 60 kN at 2/3 of [4, 8] m
assert abs(dict(r.case(1).reactions())[a][2] - 30e3) < 1e-6
assert abs(dict(r.case(2).reactions())[b][2] - 60e3 * (4 + 8 / 3) / 8) < 1e-6
