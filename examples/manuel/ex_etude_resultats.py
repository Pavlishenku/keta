# Sections d'étude imposées, déplacements, contraintes et pression du sol aux sections d'étude.
import keta

L, q = 10.0, 20e3
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
p = m.add_beam(a, b, m.add_section_rect(b=0.4, h=0.6), m.add_material(e=30e9, nu=0.2))
m.add_support(a, dx=True, dy=True, dz=True, rx=True, rz=True)
m.add_support(b, dy=True, dz=True, rz=True)
m.add_udl(p, wz=-q)
r = m.analyze(stations_abs={p: [0.0, 2.5, 5.0, 7.5, 10.0]})
cas = r.case(1)
for (xl, N, Vy, Vz, T, My, Mz), (_, u, v, w, tx, ty, tz) in zip(cas.bar_forces(p), cas.bar_deflections(p)):
    print(f"x = {xl * L:4.1f} m : My = {My / 1e3:7.2f} kN.m (q x (L-x)/2 = {q * xl * L * (L - xl * L) / 2e3:7.2f}),"
          f" w = {w * 1e3:7.3f} mm")
for xl, point, sigma, tau, vm in cas.stresses(p):
    if abs(xl - 0.5) < 1e-9:
        print(f"mi-portée, point {point:>6} : sigma = {sigma / 1e6:7.3f} MPa")

# poutre posée sur un sol élastique, sans appui vertical : w = q/k et p = q
m = keta.Model()
nd = [m.add_node(2.0 * i, 0, 0) for i in range(7)]
s = m.add_section_rect(b=1.0, h=0.5)
br = [m.add_beam(nd[i], nd[i + 1], s, m.concrete()) for i in range(6)]
m.add_support(nd[0], dx=True, dy=True, rx=True, rz=True)
for x in br:
    m.set_soil(x, kfz=5e7)
    m.add_udl(x, wz=-40e3)
c = m.analyze().case(1)
print(f"tassement = {c.displacement(nd[3])[2] * 1e3:.4f} mm (q/k = {40e3 / 5e7 * 1e3:.4f})")
print("pression du sol sous la barre 3 [kN/m] :", [round(abs(row[3]) / 1e3, 2) for row in c.soil_pressure(br[2])])
