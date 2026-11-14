# Déplacement imposé d'appui : poutre continue à trois appuis, tassement de 10 mm de l'appui central.
import keta

L, EI, d = 6.0, 30e9 * 4.0e-3, -0.010
m = keta.Model()
n = [m.add_node(L * i, 0, 0) for i in range(3)]
sec = m.add_section_props(sx=0.2, iy=4.0e-3, iz=1e-3, ix=1e-3)
b = [m.add_beam(n[i], n[i + 1], sec, m.add_material(30e9, 0.2)) for i in range(2)]
m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(n[1], dy=True, dz=True)
m.add_support(n[2], dy=True, dz=True)
m.add_support_disp(n[1], dz=d, case=1)
cas = m.analyze().case(1)
M = cas.bar_forces(b[0])[-1][5]
R = dict(cas.reactions())[n[1]][2]
print(f"déplacement imposé : {cas.displacement(n[1])[2] * 1e3:.3f} mm")
print(f"moment sur appui central : {M / 1e3:.3f} kN.m (-3 EI |d| / L² = {-3 * EI * -d / L**2 / 1e3:.3f} kN.m, fibre inférieure tendue)")
print(f"réaction de l'appui central : {R / 1e3:.3f} kN (-6 EI |d| / L³ = {-6 * EI * -d / L**3 / 1e3:.3f} kN)")
