# Non-uniform torsion and lateral-torsional buckling of an IPE 300 (warping constant Iw),
# buckling on a fixed compression state.
import math
import keta

E, NU = 210e9, 0.3
G = E / (2 * (1 + NU))
A, IY, IZ, J, IW = 53.8e-4, 8356e-8, 604e-8, 20.1e-8, 126e-9

def make_bar(L, n):
    m = keta.Model()
    st = m.add_material(e=E, nu=NU, rho=0.0)
    s = m.add_section_props(sx=A, iy=IY, iz=IZ, ix=J)
    nd = [m.add_node(L * i / n, 0, 0) for i in range(n + 1)]
    b = [m.add_beam(nd[i], nd[i + 1], s, st) for i in range(n)]
    return m, s, nd, b

# 1) Cantilever of 4 m, warping restrained at the clamp, torsional moment of 1 kN.m at the tip
L, T = 4.0, 1e3
m, s, nd, b = make_bar(L, 1)
m.fix(nd[0])
m.add_point_load(nd[-1], mx=T)
r = m.warping(case=1, Iw={s: IW}, warping_fixed=[nd[0]])
k = math.sqrt(G * J / (E * IW))
print(f"torsional rotation at the tip: {r.displacement(nd[-1])[3]:.6f} rad "
      f"(Vlasov: {T / (G * J) * (L - math.tanh(k * L) / k):.6f} rad; Saint-Venant only: {T * L / (G * J):.6f})")
forces = r.bar(b[0])
print(f"bimoment at the clamp: {abs(forces['Mw'][0]) / 1e3:.4f} kN.m2 "
      f"(T tanh(kL) / k = {T * math.tanh(k * L) / k / 1e3:.4f} kN.m2)")

# 2) Lateral-torsional buckling under uniform moment, fork supports, L = 6 m
L, M0 = 6.0, 1e3
mcr = math.pi / L * math.sqrt(E * IZ * G * J * (1 + math.pi**2 * E * IW / (G * J * L * L)))
for n in (8, 16, 32):
    m, s, nd, b = make_bar(L, n)
    m.add_support(nd[0], dx=True, dy=True, dz=True, rx=True)
    m.add_support(nd[-1], dy=True, dz=True, rx=True)
    m.add_point_load(nd[0], my=M0)
    m.add_point_load(nd[-1], my=-M0)
    lam = m.lateral_buckling(case=1, Iw={s: IW}, n_modes=1).factors[0]
    print(f"{n:3d} elements: Mcr = {lam * M0 / 1e3:.4f} kN.m (closed form {mcr / 1e3:.4f} kN.m)")

# 3) Pinned column of 5 m already compressed to 0.4 Pcr (case 1): multiplier of the variable case (case 2)
L = 5.0
pcr = math.pi**2 * E * IZ / L**2
m, s, nd, b = make_bar(L, 1)
m.add_support(nd[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(nd[-1], dy=True, dz=True)
m.add_point_load(nd[-1], fx=-0.4 * pcr, case=1)
m.add_point_load(nd[-1], fx=-1.0, case=2)
res = m.buckling_cte_var(1, 2, n_modes=1, bisection=True)
print(f"critical variable load: {res['factors'][0] / 1e3:.3f} kN (0.6 Pcr = {0.6 * pcr / 1e3:.3f} kN)")
