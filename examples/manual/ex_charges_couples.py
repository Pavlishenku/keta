# Moments on bars, partial loads and absolute abscissae.
import keta

L = 10.0
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
p = m.add_beam(a, b, m.add_section("IPE400"), m.steel())
m.add_support(a, dx=True, dy=True, dz=True, rx=True, rz=True)
m.add_support(b, dy=True, dz=True, rz=True)
m.add_bar_point_load(p, 0.5, my=50e3, frame="local", case=1)              # 50 kN.m moment at midspan
m.add_udl(p, my=3e3, frame="local", case=2)                                # distributed moment of 3 kN.m/m
m.add_udl(p, wz=-10e3, x1=0.0, x2=4.0, absolute=True, case=3)              # 10 kN/m over the first 4 metres
m.add_bar_point_load(p, 2.5, fz=-20e3, absolute=True, case=4)              # 20 kN at 2.5 m from the start
r = m.analyze()
for k, expected in ((1, 50e3 / L), (2, 3e3), (3, 10e3 * 4 * (L - 2) / L), (4, 20e3 * (L - 2.5) / L)):
    ra = dict(r.case(k).reactions())[a][2]
    print(f"case {k}: reaction at A = {ra / 1e3:8.3f} kN (by hand, absolute value: {expected / 1e3:6.3f})")
