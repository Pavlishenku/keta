# Appui biais : rouleau incliné de 30° dans le plan XZ à l'extrémité d'une poutre.
import math
import keta

m = keta.Model()
a, c = m.add_node(0, 0, 0), m.add_node(6, 0, 0)
m.add_beam(a, c, m.add_section("IPE400"), m.steel())
m.add_support(a, dx=True, dy=True, dz=True, rx=True)
t = math.radians(30)
axes = [[math.cos(t), 0.0, math.sin(t)],        # x local : direction libre du rouleau
        [0.0, 1.0, 0.0],                          # y local
        [-math.sin(t), 0.0, math.cos(t)]]         # z local : normale au plan de roulement
m.add_skew_support(c, axes, dy=True, dz=True)     # bloque y et z locaux
m.add_point_load(c, fz=-50e3)

cas = m.analyze().case(1)
noeud, loc, glob = cas.skew_reactions()[0]
print(f"noeud {noeud} : réaction locale  Rz' = {loc[2] / 1e3:.3f} kN")
print(f"           réaction globale Rx = {glob[0] / 1e3:.3f} kN, Rz = {glob[2] / 1e3:.3f} kN")
print(f"réaction de l'appui fixe : Rx = {cas.reactions()[0][1][0] / 1e3:.3f} kN")
assert abs(loc[2] - 50e3 / math.cos(t)) < 1e-6 * 50e3
