# Orientation de la section (angle beta) : une même console IPE 300 fléchie
# autour de son axe fort (beta = 0) puis de son axe faible (beta = pi/2).
import math
import keta

def z(x):                     # supprime le zéro négatif d'affichage
    return 0.0 if abs(x) < 1e-6 else x

for beta in (0.0, math.pi / 2):
    m = keta.Model()
    a, b = m.add_node(0, 0, 0), m.add_node(3, 0, 0)
    m.add_beam(a, b, m.add_section("IPE300"), m.steel(), beta=beta)
    m.fix(a)
    m.add_point_load(b, fz=-5e3)
    cas = m.analyze().case(1)
    uz = cas.displacement(b)[2]
    x, N, Vy, Vz, T, My, Mz = cas.bar_forces(1)[0]
    print(f"beta = {beta:.4f} rad : uz = {uz * 1e3:8.3f} mm, My(0) = {z(My) / 1e3:7.3f} kN.m, Mz(0) = {z(Mz) / 1e3:7.3f} kN.m")

# Repère local d'une barre verticale et d'une barre inclinée (lecture par la flèche d'une charge)
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(0, 0, 4)
m.add_beam(a, b, m.add_section_rect(0.3, 0.5), m.concrete())
m.fix(a)
m.add_point_load(b, fx=10e3)          # force horizontale en tête de poteau
x, N, Vy, Vz, T, My, Mz = m.analyze().case(1).bar_forces(1)[0]
print(f"poteau, pied : Vy = {z(Vy) / 1e3:.3f} kN, Vz = {Vz / 1e3:.3f} kN, My = {My / 1e3:.3f} kN.m, Mz = {z(Mz) / 1e3:.3f} kN.m")
