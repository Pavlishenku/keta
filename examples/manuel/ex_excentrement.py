# Excentrements : tablier mixte schématique, poutre décalée de 0,6 m sous la ligne des noeuds.
import keta

def z(x):                     # supprime le zéro négatif d'affichage
    return 0.0 if abs(x) < 1e-6 else x

L = 10.0
for e in (0.0, -0.6):
    m = keta.Model()
    a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
    poutre = m.add_beam(a, b, m.add_section("HEB500"), m.steel(), ecc=(0.0, 0.0, e))
    m.add_support(a, dx=True, dy=True, dz=True, rx=True)
    m.add_support(b, dy=True, dz=True)
    m.add_point_load(b, fx=-200e3)                  # poussée horizontale sur l'appui mobile
    cas = m.analyze().case(1)
    st = cas.bar_forces(poutre)
    print(f"e = {e:5.2f} m : N = {st[5][1] / 1e3:8.2f} kN, My(milieu) = {z(st[5][5]) / 1e3:8.2f} kN.m, "
          f"ux(B) = {cas.displacement(b)[0] * 1e3:7.4f} mm, ry(B) = {z(cas.displacement(b)[4]) * 1e3:7.4f} mrad")
