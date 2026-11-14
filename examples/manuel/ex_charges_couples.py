# Couples sur barres, charges partielles et abscisses absolues.
import keta

L = 10.0
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
p = m.add_beam(a, b, m.add_section("IPE400"), m.steel())
m.add_support(a, dx=True, dy=True, dz=True, rx=True, rz=True)
m.add_support(b, dy=True, dz=True, rz=True)
m.add_bar_point_load(p, 0.5, my=50e3, frame="local", case=1)              # couple de 50 kN.m à mi-portée
m.add_udl(p, my=3e3, frame="local", case=2)                                # couple réparti de 3 kN.m/m
m.add_udl(p, wz=-10e3, x1=0.0, x2=4.0, absolute=True, case=3)              # 10 kN/m sur les 4 premiers mètres
m.add_bar_point_load(p, 2.5, fz=-20e3, absolute=True, case=4)              # 20 kN à 2,5 m de l'origine
r = m.analyze()
for k, attendu in ((1, 50e3 / L), (2, 3e3), (3, 10e3 * 4 * (L - 2) / L), (4, 20e3 * (L - 2.5) / L)):
    ra = dict(r.case(k).reactions())[a][2]
    print(f"cas {k} : réaction en A = {ra / 1e3:8.3f} kN (main, en valeur absolue : {attendu / 1e3:6.3f})")
