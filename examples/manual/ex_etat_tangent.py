# Saved tangent state: pinned-pinned column of 10 m (EI = 210e9 x 8e-6) compressed to 30 % of the
# Euler load. Modes about the loaded state, then buckling of an additional load.
import math
import keta

L, E, I, A, RHO = 10.0, 210e9, 8e-6, 4e-3, 7850.0
PE = math.pi ** 2 * E * I / L ** 2                      # Euler load
m = keta.Model()
mat = m.add_material(e=E, nu=0.3, rho=RHO)
sec = m.add_section_props(sx=A, iy=I, iz=50 * I, ix=1e-5)
nd = [m.add_node(L * i / 20, 0, 0) for i in range(21)]
for i in range(20):
    m.add_beam(nd[i], nd[i + 1], sec, mat)
m.add_support(nd[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(nd[20], dy=True, dz=True)
m.add_point_load(nd[20], fx=-0.3 * PE, case=1)          # state: compression 0.3 PE

t = m.tangent(case=1)                                     # tangent K + Kg(N) of the state
print(f"N in bar 1: {t.axial(1) / 1e3:.2f} kN (PE = {PE / 1e3:.2f} kN)")
w0 = m.modal(n_modes=1).omegas[0]
w = t.modal(n_modes=1).omegas[0]
print(f"omega1 unloaded {w0:.4f} rad/s, under the state {w:.4f} rad/s, ratio of the squares {w * w / w0 ** 2:.4f} "
      f"(1 - P/PE = 0.7000)")
m.add_point_load(nd[20], fx=-1e3, case=2)               # additional load of 1 kN
b = t.buckling(m, case=2)
print(f"multiplier of the additional load: {b.factors[0]:.2f} ((PE - P0)/P1 = {0.7 * PE / 1e3:.2f})")
