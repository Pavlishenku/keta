# Special staged analyses: construction in large displacements (segment placed on the
# tangent of the deformed cantilever) and fibre column creeping between two stages.
import math
import keta

# 1) Cantilever of 4 m built in two segments of 2 m; a moment at the end of the first segment
m = keta.Model()
mat = m.add_material(e=2e11, nu=0.3)
sec = m.add_section_props(sx=1e-2, iy=1e-5, iz=1e-5, ix=2e-5)
nd = [m.add_node(0.5 * i, 0, 0) for i in range(9)]
b = [m.add_beam(nd[i], nd[i + 1], sec, mat) for i in range(8)]
m.fix(nd[0])
for x in nd[1:]:
    m.add_support(x, dy=True, rx=True, rz=True)
m.add_point_load(nd[4], my=0.2 * 2e6 / 2.0, case=1)        # end rotation of 0.2 rad
states = m.staged_large_displacements([dict(activate=b[:4], case=1), dict(activate=b[4:])])
print(f"stage 1: rotation at x = 2 m: {states[0][nd[4]][4]:.6f} rad, deflection {states[0][nd[4]][2] * 1e3:.3f} mm")
print(f"stage 2: rotation at x = 4 m: {states[1][nd[8]][4]:.6f} rad, deflection {states[1][nd[8]][2] * 1e3:.3f} mm")

# 2) Fibre column of 3 m (section 0.30 x 0.30 m), loaded to 9 MPa, Kelvin creep
#    (a = 2/E, tau = 60 days) for 300 days: eps = s/E + s a (1 - exp(-t / tau))
m = keta.Model()
c = m.add_material(e=3e10, nu=0.2)
s = m.add_section_props(sx=0.09, iy=0.09**2 / 12, iz=0.09**2 / 12, ix=1e-3)
a_, b_ = m.add_node(0, 0, 0), m.add_node(0, 0, 3)
column = m.add_beam(a_, b_, s, c)
m.fix(a_)
m.add_support(b_, dx=True, dy=True, rx=True, ry=True, rz=True)
m.add_point_load(b_, fz=-0.09 * 9e6)
law = dict(kind="elastic", E=3e10, creep=dict(a=[2 / 3e10], tau=[60.0]))
r = m.fiber_frame({column: dict(rect=(0.3, 0.3, 2, 2), law=law)},
                  stages=[dict(activate=[column], case=1, duration=300.0, time_steps=6)])
eps = r["node_disp"][b_][2] / 3
print(f"fibres: strain at 300 days = {eps * 1e6:.3f} µm/m; "
      f"closed form {-9e6 / 3e10 * (1 + 2 * (1 - math.exp(-300 / 60))) * 1e6:.3f} µm/m")
print("dates of the trace [days]:", [round(t, 1) for t in r["time"]])
