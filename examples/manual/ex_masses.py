# Nodal masses and distributed masses added to bars.
import math
import keta

E, I, L, M = 210e9, 1e-4, 4.0, 500.0
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
sec = m.add_section_props(sx=0.01, iy=I, iz=I, ix=I)
m.add_beam(a, b, sec, m.add_material(e=E, nu=0.3))          # material without density
m.fix(a)
m.add_nodal_mass(b, diag=[M, M, M, 0, 0, 0])
f = m.modal(n_modes=1).frequencies[0]
print(f"tip mass: f = {f:.4f} Hz (sqrt(3EI/L3/M)/2pi = {math.sqrt(3 * E * I / L**3 / M) / (2 * math.pi):.4f})")

# simply supported beam of 20 m, massless material, added distributed mass of 3 t/m
L, Ib, mu = 20.0, 0.5, 3000.0
m = keta.Model()
nd = [m.add_node(L * i / 40, 0, 0) for i in range(41)]
sec = m.add_section_props(sx=1.5, iy=Ib, iz=2.0, ix=0.3)
mat = m.add_material(e=30e9, nu=0.2)
bars = [m.add_beam(nd[i], nd[i + 1], sec, mat) for i in range(40)]
m.add_support(nd[0], dx=True, dy=True, dz=True, rx=True, rz=True)
m.add_support(nd[-1], dy=True, dz=True, rz=True)
for x in nd[1:-1]:
    m.add_support(x, dy=True, rx=True, rz=True)
for br in bars:
    m.add_bar_mass(br, "uniform", m1=mu)
f1 = m.modal(n_modes=1).frequencies[0]
print(f"distributed mass: f1 = {f1:.4f} Hz (pi/2L2 sqrt(EI/mu) = {math.pi / (2 * L * L) * math.sqrt(30e9 * Ib / mu):.4f})")
