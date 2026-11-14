# Influence lines: continuous bridge with three spans 25 + 30 + 25 m, moment over a pier
# and moment at the middle of the central span; extremes under point load, distributed load and vehicles.
import keta

spans = [25.0, 30.0, 25.0]
m = keta.Model()
xs = [0.0]
for L in spans:
    xs += [xs[-1] + L * i / 10 for i in range(1, 11)]
n = [m.add_node(x, 0, 0) for x in xs]
sec = m.add_section_props(sx=6.0, iy=4.0, iz=40.0, ix=5.0)
b = [m.add_beam(n[i], n[i + 1], sec, m.concrete()) for i in range(len(n) - 1)]
m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
for i in (10, 20, 30):
    m.add_support(n[i], dy=True, dz=True)

# My at the end of bar 10 (over pier 1) and at the middle of the central span (end of bar 15)
line_pier = m.influence_line(b, dict(kind="bar_force", bar=b[9], x=2.5, comp=4))
line_mid = m.influence_line(b, dict(kind="bar_force", bar=b[14], x=3.0, comp=4))
print(f"line: {len(line_pier.s)} abscissae from {line_pier.s[0]:.1f} to {line_pier.s[-1]:.1f} m")
for name, line in (("M over pier", line_pier), ("M midspan", line_mid)):
    mx, mn = line.point_extreme()
    print(f"{name:12}: unit load -> max {mx:7.3f}, min {mn:7.3f} m; "
          f"areas + {line.integral_pos():8.3f} / - {line.integral_neg():8.3f} m2")
    print(f"{'':12}  positive zones: " + ", ".join(f"[{a:.2f}; {c:.2f}]" for a, c, _ in line.zones(1)))
    print(f"{'':12}  ordinate at s = 40 m: {line.at(40.0):.4f} m")

# Moment at the middle of the central span: vehicles and distributed loads (characteristic values, N.m)
q = 9e3 * 3.0                                    # 9 kN/m2 on a lane of 3 m
cases = [("LM1 tandem lane 1", keta.convoy_extreme(line_mid, keta.lm1_tandem(1), 2000)),
         ("distributed 27 kN/m", keta.udl_extreme(line_mid, q)),
         ("LM71", keta.convoy_extreme(line_mid, keta.rail_lm71(), 2000))]
for name, (mx, mn, rmx, rmn) in cases:
    print(f"{name:19}: max = {mx / 1e3:9.1f} kN.m, min = {mn / 1e3:9.1f} kN.m")
