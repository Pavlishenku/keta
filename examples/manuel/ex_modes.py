# Modes propres : poutre isostatique de 20 m (masse répartie), puis ajout d'une masse ponctuelle.
import math
import keta

L, E, I, rho, A = 20.0, 35e9, 0.5, 2500.0, 1.5

def poutre(masse=0.0):
    m = keta.Model()
    nd = [m.add_node(L * i / 20, 0, 0) for i in range(21)]
    sec = m.add_section_props(sx=A, iy=I, iz=2.0, ix=0.3)
    mat = m.add_material(e=E, nu=0.2, rho=rho)
    for i in range(20):
        m.add_beam(nd[i], nd[i + 1], sec, mat)
    m.add_support(nd[0], dx=True, dy=True, dz=True, rx=True, rz=True)
    m.add_support(nd[-1], dy=True, dz=True, rz=True)
    for x in nd[1:-1]:
        m.add_support(x, dx=True, dy=True, rx=True, rz=True)   # flexion dans le plan XZ seulement
    if masse:
        m.add_point_mass(nd[10], masse)
    return m, nd

m, nd = poutre()
mo = m.modal(n_modes=4)
print(mo)
mu = rho * A
for k, f in enumerate(mo.frequencies, start=1):
    fth = k * k * math.pi / (2 * L * L) * math.sqrt(E * I / mu)
    print(f"mode {k} : f = {f:8.4f} Hz  T = {mo.periods[k - 1]:.4f} s  (Euler-Bernoulli : {fth:8.4f} Hz)")
print("masse totale par direction [t] :", [round(x / 1e3, 3) for x in mo.total_mass()[:3]])
for k, r in enumerate(mo.mass_ratio(), start=1):
    print(f"mode {k} : masse effective selon Z = {mo.effective_mass()[k - 1][2] / 1e3:8.3f} t ({100 * r[2]:6.2f} %)")
x, phi = mo.mode(0)[10]
print(f"déformée modale 1 au noeud {x} (amplitude unité en crête) : uz = {phi[2]:.4f}")

m2, nd2 = poutre(masse=20e3)
f2 = m2.modal(n_modes=1, method="dense").frequencies[0]
print(f"avec une masse de 20 t à mi-portée : f1 = {f2:.4f} Hz")
