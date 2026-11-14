# Phasages particuliers : construction en grands déplacements (segment posé sur la
# tangente de la console déformée) et poteau à fibres fluant entre deux phases.
import math
import keta

# 1) Console de 4 m construite en deux segments de 2 m ; un couple en bout du premier segment
m = keta.Model()
mat = m.add_material(e=2e11, nu=0.3)
sec = m.add_section_props(sx=1e-2, iy=1e-5, iz=1e-5, ix=2e-5)
nd = [m.add_node(0.5 * i, 0, 0) for i in range(9)]
b = [m.add_beam(nd[i], nd[i + 1], sec, mat) for i in range(8)]
m.fix(nd[0])
for x in nd[1:]:
    m.add_support(x, dy=True, rx=True, rz=True)
m.add_point_load(nd[4], my=0.2 * 2e6 / 2.0, case=1)        # rotation d'extrémité 0,2 rad
etats = m.staged_large_displacements([dict(activate=b[:4], case=1), dict(activate=b[4:])])
print(f"phase 1 : rotation en x = 2 m : {etats[0][nd[4]][4]:.6f} rad, flèche {etats[0][nd[4]][2] * 1e3:.3f} mm")
print(f"phase 2 : rotation en x = 4 m : {etats[1][nd[8]][4]:.6f} rad, flèche {etats[1][nd[8]][2] * 1e3:.3f} mm")

# 2) Poteau de 3 m en fibres (section 0,30 x 0,30 m), chargé à 9 MPa, fluage de Kelvin
#    (a = 2/E, tau = 60 j) pendant 300 jours : eps = s/E + s a (1 - exp(-t / tau))
m = keta.Model()
c = m.add_material(e=3e10, nu=0.2)
s = m.add_section_props(sx=0.09, iy=0.09**2 / 12, iz=0.09**2 / 12, ix=1e-3)
a_, b_ = m.add_node(0, 0, 0), m.add_node(0, 0, 3)
poteau = m.add_beam(a_, b_, s, c)
m.fix(a_)
m.add_support(b_, dx=True, dy=True, rx=True, ry=True, rz=True)
m.add_point_load(b_, fz=-0.09 * 9e6)
loi = dict(kind="elastic", E=3e10, creep=dict(a=[2 / 3e10], tau=[60.0]))
r = m.fiber_frame({poteau: dict(rect=(0.3, 0.3, 2, 2), law=loi)},
                  stages=[dict(activate=[poteau], case=1, duration=300.0, time_steps=6)])
eps = r["node_disp"][b_][2] / 3
print(f"fibres : déformation à 300 j = {eps * 1e6:.3f} µm/m ; "
      f"forme fermée {-9e6 / 3e10 * (1 + 2 * (1 - math.exp(-300 / 60))) * 1e6:.3f} µm/m")
print("dates de la trace [j] :", [round(t, 1) for t in r["time"]])
