# Sections quelconques : polygone avec évidement (éléments finis de section) et parois minces.
import math
import keta

# rectangle plein 0,8 x 0,4 m : la torsion comparée à la série de Saint-Venant
b, h = 0.8, 0.4
rect = [(-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2)]
p = keta.section_polygon(rect, mesh_size=0.02)
serie = sum(math.tanh(k * math.pi * b / (2 * h)) / k**5 for k in range(1, 200, 2))
J = b * h**3 / 3 * (1 - 192 / math.pi**5 * h / b * serie)
print(f"A = {p['A']:.4f} m2, Iy = {p['Iy']:.6f} m4 (b h3/12 = {b * h**3 / 12:.6f})")
print(f"J = {p['J']:.6f} m4 (série : {J:.6f}), Avz = {p['Avz']:.4f} m2")

# caisson évidé : le noyau calcule J, Iw et le centre de cisaillement
vide = [(-0.3, -0.1), (0.3, -0.1), (0.3, 0.1), (-0.3, 0.1)]
q = keta.section_polygon(rect, [vide], mesh_size=0.02)
print(f"caisson : A = {q['A']:.4f} m2, Iy = {q['Iy']:.6f} m4, J = {q['J']:.6f} m4")

# parois minces : caisson carré de 1 m, épaisseur 10 mm (Bredt : J = a3 t)
a, t = 1.0, 0.01
c = keta.thinwall_section([(0, 0), (a, 0), (a, a), (0, a)], [(0, 1, t), (1, 2, t), (2, 3, t), (3, 0, t)])
print(f"Bredt : J = {c['J']:.6f} m4 (a3 t = {a**3 * t:.6f}), {c['cells']} cellule")

# section polygonale dans un modèle : console en torsion
m = keta.Model()
s = m.add_section_polygon(rect, mesh_size=0.02)
n1, n2 = m.add_node(0, 0, 0), m.add_node(4, 0, 0)
m.add_beam(n1, n2, s, m.add_material(e=30e9, nu=0.2))
m.fix(n1)
m.add_point_load(n2, mx=1e3)
rx = m.analyze().case(1).displacement(n2)[3]
print(f"rotation de torsion = {rx * 1e6:.4f} microrad (T L / G J = {1e3 * 4 / (12.5e9 * J) * 1e6:.4f})")
