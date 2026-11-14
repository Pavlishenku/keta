# General sections: polygon with a void (section finite elements) and thin walls.
import math
import keta

# solid rectangle 0.8 x 0.4 m: torsion compared with the Saint-Venant series
b, h = 0.8, 0.4
rect = [(-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2)]
p = keta.section_polygon(rect, mesh_size=0.02)
series = sum(math.tanh(k * math.pi * b / (2 * h)) / k**5 for k in range(1, 200, 2))
J = b * h**3 / 3 * (1 - 192 / math.pi**5 * h / b * series)
print(f"A = {p['A']:.4f} m2, Iy = {p['Iy']:.6f} m4 (b h3/12 = {b * h**3 / 12:.6f})")
print(f"J = {p['J']:.6f} m4 (series: {J:.6f}), Avz = {p['Avz']:.4f} m2")

# hollow box: the kernel computes J, Iw and the shear centre
void = [(-0.3, -0.1), (0.3, -0.1), (0.3, 0.1), (-0.3, 0.1)]
q = keta.section_polygon(rect, [void], mesh_size=0.02)
print(f"box: A = {q['A']:.4f} m2, Iy = {q['Iy']:.6f} m4, J = {q['J']:.6f} m4")

# thin walls: square box of 1 m, thickness 10 mm (Bredt: J = a3 t)
a, t = 1.0, 0.01
c = keta.thinwall_section([(0, 0), (a, 0), (a, a), (0, a)], [(0, 1, t), (1, 2, t), (2, 3, t), (3, 0, t)])
print(f"Bredt: J = {c['J']:.6f} m4 (a3 t = {a**3 * t:.6f}), {c['cells']} cell")

# polygonal section in a model: cantilever in torsion
m = keta.Model()
s = m.add_section_polygon(rect, mesh_size=0.02)
n1, n2 = m.add_node(0, 0, 0), m.add_node(4, 0, 0)
m.add_beam(n1, n2, s, m.add_material(e=30e9, nu=0.2))
m.fix(n1)
m.add_point_load(n2, mx=1e3)
rx = m.analyze().case(1).displacement(n2)[3]
print(f"torsional rotation = {rx * 1e6:.4f} microrad (T L / G J = {1e3 * 4 / (12.5e9 * J) * 1e6:.4f})")
