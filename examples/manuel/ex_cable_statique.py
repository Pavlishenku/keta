# Câble de précontrainte dans un cas de charge statique : tracé parabolique exact, poutre isostatique.
import keta

L, P, f = 20.0, 2e6, 0.5
m = keta.Model()
nd = [m.add_node(L * i / 8, 0, 0) for i in range(9)]
s = m.add_section_rect(b=0.6, h=1.2)
b = [m.add_beam(nd[i], nd[i + 1], s, m.concrete()) for i in range(8)]
m.add_support(nd[0], dx=True, dy=True, dz=True, rx=True, rz=True)
m.add_support(nd[-1], dy=True, dz=True, rz=True)
# e(s) = -4 f s (L - s) / L2 : un seul segment cubique (c0, c1, c2, c3)
cable = m.add_tendon(b, P0=P, segments=[dict(s0=0.0, s1=L, cz=(0.0, -4 * f / L, 4 * f / L**2))],
                     mu=0.0, k=0.0, slip=0.0)
print("charges équivalentes ajoutées :", m.add_tendon_loads(cable, case=1))
c = m.analyze().case(1)
mi = c.bar_forces(b[3])[-1]
print(f"mi-portée : N = {mi[1] / 1e3:.1f} kN, My = {mi[5] / 1e3:.1f} kN.m (P f = {P * f / 1e3:.1f})")
print(f"réaction verticale en A = {abs(dict(c.reactions())[nd[0]][2]):.3f} N (isostatique : 0)")
