# Câbles en chaînette élastique : câble isolé (comparaison à la forme fermée d'Irvine),
# puis mât haubané résolu en statique couplée câbles + ossature et en grands déplacements.
import math
import keta

# 1) Câble de 100 m de portée horizontale, poids 50 N/m, tension horizontale imposée 20 kN
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(100, 0, 0)
m.fix(a)
m.fix(b)
c = m.add_cable(a, b, ea=2.0e8, w=50.0, h=20e3)
m.add_point_load(b, fz=0.0)                     # un cas de charge (vide) est nécessaire
cs = m.cable_static()
H, V1, V2, Tmax, L0 = cs.cable(c)
forme = cs.cable_shape(c)
fl = min(p[2] for p in forme)
print(f"H = {H / 1e3:.4f} kN, V1 = {V1 / 1e3:.4f} kN, Tmax = {Tmax / 1e3:.4f} kN, L0 = {L0:.4f} m")
print(f"flèche (21 points) {fl:.4f} m ; chaînette inextensible H/w (cosh(wL/2H) - 1) = "
      f"{-H / 50 * (math.cosh(50 * 100 / (2 * H)) - 1):.4f} m")

# 2) Mât de 6 m encastré, console de 8 m, hauban tendu (H = 20 kN) de la tête du mât au bout de console
m = keta.Model()
pied, tete, bout = m.add_node(0, 0, 0), m.add_node(0, 0, 6), m.add_node(8, 0, 0)
acier, heb = m.steel(), m.add_section("HEB200")
mat = m.add_beam(pied, tete, heb, acier)
console = m.add_beam(pied, bout, heb, acier)
m.fix(pied)
h = m.add_cable(tete, bout, ea=2.0e8, w=50.0, h=2.0e4)
m.add_point_load(bout, fz=-5000.0)
cs = m.cable_static(substeps=2)
print(cs)
H, V1, V2, Tmax, L0 = cs.cable(h)
print(f"statique câbles : H = {H / 1e3:.3f} kN, Tmax = {Tmax / 1e3:.3f} kN, L0 = {L0:.5f} m, "
      f"flèche du bout {cs.displacement(bout)[2] * 1e3:.3f} mm")
print(f"moment en pied de console : {cs.bar_forces(console)[0][5] / 1e3:.3f} kN.m")
nl = m.nonlinear(substeps=2)
t1, t2 = nl.cable_tension(h)
print(f"grands déplacements : tensions d'extrémité {t1 / 1e3:.3f} / {t2 / 1e3:.3f} kN, "
      f"flèche du bout {nl.displacement(bout)[2] * 1e3:.3f} mm")
