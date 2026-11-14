# Offsets specific to each end, translation and torsion releases.
import keta

def z(x):
    return 0.0 if abs(x) < 1e-9 else x

# 1) 6 m cantilever whose axis is raised by 0.4 m at both ends: H in the line of the nodes
L, e, H = 6.0, 0.4, 100e3
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
sec = m.add_section_props(sx=0.3, iy=0.01, iz=0.05, ix=0.01)
mat = m.add_material(e=30e9, nu=0.2)
beam = m.add_beam(a, b, sec, mat)
m.set_eccentricity(beam, origin=(0, 0, e), end=(0, 0, e))
m.fix(a)
m.add_point_load(b, fx=H)
case_ = m.analyze().case(1)
ux = case_.displacement(b)[0]
print(f"ux(B) = {ux * 1e3:.4f} mm  (H L/EA + H e2 L/EI = {(H * L / (30e9 * 0.3) + H * e * e * L / (30e9 * 0.01)) * 1e3:.4f} mm)")
print(f"My in the bar = {case_.bar_forces(beam)[5][5] / 1e3:.2f} kN.m  (H e = {H * e / 1e3:.2f})")

# 2) same cantilever, only the end is raised: the bar is inclined
m.set_eccentricity(beam, end=(0, 0, e))
case_ = m.analyze().case(1)
print(f"clamping moment = {z(case_.reactions()[0][1][4]):.2f} N.m, "
      f"moment at the bar end = {case_.bar_forces(beam)[-1][5] / 1e3:.2f} kN.m")

# 3) two aligned bars clamped at both ends; the first one released in dx and rx at its end
m = keta.Model()
n = [m.add_node(5.0 * i, 0, 0) for i in range(3)]
s = m.add_section_rect(b=0.3, h=0.5)
b1 = m.add_beam(n[0], n[1], s, m.concrete())
m.add_beam(n[1], n[2], s, m.concrete())
m.fix(n[0]); m.fix(n[2])
m.set_releases(b1, end="dx rx")
m.add_point_load(n[1], fx=10e3, mx=2e3)
for nd, r in m.analyze().case(1).reactions():
    print(f"node {nd}: Fx = {z(r[0]) / 1e3:7.2f} kN, Mx = {z(r[3]) / 1e3:6.2f} kN.m")
