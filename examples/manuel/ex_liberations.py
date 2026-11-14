# Excentrements propres à chaque extrémité, libérations de translation et de torsion.
import keta

def z(x):
    return 0.0 if abs(x) < 1e-9 else x

# 1) console de 6 m dont l'axe est relevé de 0,4 m aux deux extrémités : H dans la ligne des noeuds
L, e, H = 6.0, 0.4, 100e3
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
sec = m.add_section_props(sx=0.3, iy=0.01, iz=0.05, ix=0.01)
mat = m.add_material(e=30e9, nu=0.2)
poutre = m.add_beam(a, b, sec, mat)
m.set_eccentricity(poutre, origin=(0, 0, e), end=(0, 0, e))
m.fix(a)
m.add_point_load(b, fx=H)
cas = m.analyze().case(1)
ux = cas.displacement(b)[0]
print(f"ux(B) = {ux * 1e3:.4f} mm  (H L/EA + H e2 L/EI = {(H * L / (30e9 * 0.3) + H * e * e * L / (30e9 * 0.01)) * 1e3:.4f} mm)")
print(f"My dans la barre = {cas.bar_forces(poutre)[5][5] / 1e3:.2f} kN.m  (H e = {H * e / 1e3:.2f})")

# 2) même console, seule l'extrémité est relevée : la barre est inclinée
m.set_eccentricity(poutre, end=(0, 0, e))
cas = m.analyze().case(1)
print(f"moment d'encastrement = {z(cas.reactions()[0][1][4]):.2f} N.m, "
      f"moment en bout de barre = {cas.bar_forces(poutre)[-1][5] / 1e3:.2f} kN.m")

# 3) deux barres alignées encastrées aux deux bouts ; la première libérée en dx et rx à son extrémité
m = keta.Model()
n = [m.add_node(5.0 * i, 0, 0) for i in range(3)]
s = m.add_section_rect(b=0.3, h=0.5)
b1 = m.add_beam(n[0], n[1], s, m.concrete())
m.add_beam(n[1], n[2], s, m.concrete())
m.fix(n[0]); m.fix(n[2])
m.set_releases(b1, end="dx rx")
m.add_point_load(n[1], fx=10e3, mx=2e3)
for nd, r in m.analyze().case(1).reactions():
    print(f"noeud {nd} : Fx = {z(r[0]) / 1e3:7.2f} kN, Mx = {z(r[3]) / 1e3:6.2f} kN.m")
