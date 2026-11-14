# Mode shapes: simply supported beam of 20 m (distributed mass), then addition of a point mass.
import math
import keta

L, E, I, rho, A = 20.0, 35e9, 0.5, 2500.0, 1.5

def make_beam(mass=0.0):
    m = keta.Model()
    nd = [m.add_node(L * i / 20, 0, 0) for i in range(21)]
    sec = m.add_section_props(sx=A, iy=I, iz=2.0, ix=0.3)
    mat = m.add_material(e=E, nu=0.2, rho=rho)
    for i in range(20):
        m.add_beam(nd[i], nd[i + 1], sec, mat)
    m.add_support(nd[0], dx=True, dy=True, dz=True, rx=True, rz=True)
    m.add_support(nd[-1], dy=True, dz=True, rz=True)
    for x in nd[1:-1]:
        m.add_support(x, dx=True, dy=True, rx=True, rz=True)   # bending in the XZ plane only
    if mass:
        m.add_point_mass(nd[10], mass)
    return m, nd

m, nd = make_beam()
mo = m.modal(n_modes=4)
print(mo)
mu = rho * A
for k, f in enumerate(mo.frequencies, start=1):
    fth = k * k * math.pi / (2 * L * L) * math.sqrt(E * I / mu)
    print(f"mode {k}: f = {f:8.4f} Hz  T = {mo.periods[k - 1]:.4f} s  (Euler-Bernoulli: {fth:8.4f} Hz)")
print("total mass per direction [t]:", [round(x / 1e3, 3) for x in mo.total_mass()[:3]])
for k, r in enumerate(mo.mass_ratio(), start=1):
    print(f"mode {k}: effective mass along Z = {mo.effective_mass()[k - 1][2] / 1e3:8.3f} t ({100 * r[2]:6.2f} %)")
x, phi = mo.mode(0)[10]
print(f"mode shape 1 at node {x} (unit peak amplitude): uz = {phi[2]:.4f}")

m2, nd2 = make_beam(mass=20e3)
f2 = m2.modal(n_modes=1, method="dense").frequencies[0]
print(f"with a mass of 20 t at midspan: f1 = {f2:.4f} Hz")
