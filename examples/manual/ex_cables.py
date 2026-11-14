# Elastic catenary cables: single cable (compared with Irvine's closed form),
# then a guyed mast solved by coupled cable + frame statics and in large displacements.
import math
import keta

# 1) Cable of 100 m horizontal span, weight 50 N/m, imposed horizontal tension 20 kN
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(100, 0, 0)
m.fix(a)
m.fix(b)
c = m.add_cable(a, b, ea=2.0e8, w=50.0, h=20e3)
m.add_point_load(b, fz=0.0)                     # a (empty) load case is required
cs = m.cable_static()
H, V1, V2, Tmax, L0 = cs.cable(c)
shape = cs.cable_shape(c)
sag = min(p[2] for p in shape)
print(f"H = {H / 1e3:.4f} kN, V1 = {V1 / 1e3:.4f} kN, Tmax = {Tmax / 1e3:.4f} kN, L0 = {L0:.4f} m")
print(f"sag (21 points) {sag:.4f} m; inextensible catenary H/w (cosh(wL/2H) - 1) = "
      f"{-H / 50 * (math.cosh(50 * 100 / (2 * H)) - 1):.4f} m")

# 2) Clamped 6 m mast, 8 m cantilever, stay tensioned (H = 20 kN) from the mast head to the cantilever tip
m = keta.Model()
foot, head, tip = m.add_node(0, 0, 0), m.add_node(0, 0, 6), m.add_node(8, 0, 0)
steel, heb = m.steel(), m.add_section("HEB200")
mast = m.add_beam(foot, head, heb, steel)
cantilever = m.add_beam(foot, tip, heb, steel)
m.fix(foot)
h = m.add_cable(head, tip, ea=2.0e8, w=50.0, h=2.0e4)
m.add_point_load(tip, fz=-5000.0)
cs = m.cable_static(substeps=2)
print(cs)
H, V1, V2, Tmax, L0 = cs.cable(h)
print(f"cable statics: H = {H / 1e3:.3f} kN, Tmax = {Tmax / 1e3:.3f} kN, L0 = {L0:.5f} m, "
      f"tip deflection {cs.displacement(tip)[2] * 1e3:.3f} mm")
print(f"moment at the cantilever root: {cs.bar_forces(cantilever)[0][5] / 1e3:.3f} kN.m")
nl = m.nonlinear(substeps=2)
t1, t2 = nl.cable_tension(h)
print(f"large displacements: end tensions {t1 / 1e3:.3f} / {t2 / 1e3:.3f} kN, "
      f"tip deflection {nl.displacement(tip)[2] * 1e3:.3f} mm")
