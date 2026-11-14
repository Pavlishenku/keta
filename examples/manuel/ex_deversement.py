# Torsion non uniforme et déversement d'un IPE 300 (inertie de gauchissement Iw),
# flambement sur un état de compression fixé.
import math
import keta

E, NU = 210e9, 0.3
G = E / (2 * (1 + NU))
A, IY, IZ, J, IW = 53.8e-4, 8356e-8, 604e-8, 20.1e-8, 126e-9

def barre(L, n):
    m = keta.Model()
    st = m.add_material(e=E, nu=NU, rho=0.0)
    s = m.add_section_props(sx=A, iy=IY, iz=IZ, ix=J)
    nd = [m.add_node(L * i / n, 0, 0) for i in range(n + 1)]
    b = [m.add_beam(nd[i], nd[i + 1], s, st) for i in range(n)]
    return m, s, nd, b

# 1) Console de 4 m, gauchissement empêché à l'encastrement, couple de torsion 1 kN.m en bout
L, T = 4.0, 1e3
m, s, nd, b = barre(L, 1)
m.fix(nd[0])
m.add_point_load(nd[-1], mx=T)
r = m.warping(case=1, Iw={s: IW}, warping_fixed=[nd[0]])
k = math.sqrt(G * J / (E * IW))
print(f"rotation de torsion en bout : {r.displacement(nd[-1])[3]:.6f} rad "
      f"(Vlasov : {T / (G * J) * (L - math.tanh(k * L) / k):.6f} rad ; Saint-Venant seul : {T * L / (G * J):.6f})")
efforts = r.bar(b[0])
print(f"bimoment à l'encastrement : {abs(efforts['Mw'][0]) / 1e3:.4f} kN.m2 "
      f"(T tanh(kL) / k = {T * math.tanh(k * L) / k / 1e3:.4f} kN.m2)")

# 2) Déversement sous moment uniforme, appuis à fourche, L = 6 m
L, M0 = 6.0, 1e3
mcr = math.pi / L * math.sqrt(E * IZ * G * J * (1 + math.pi**2 * E * IW / (G * J * L * L)))
for n in (8, 16, 32):
    m, s, nd, b = barre(L, n)
    m.add_support(nd[0], dx=True, dy=True, dz=True, rx=True)
    m.add_support(nd[-1], dy=True, dz=True, rx=True)
    m.add_point_load(nd[0], my=M0)
    m.add_point_load(nd[-1], my=-M0)
    lam = m.lateral_buckling(case=1, Iw={s: IW}, n_modes=1).factors[0]
    print(f"{n:3d} éléments : Mcr = {lam * M0 / 1e3:.4f} kN.m (forme fermée {mcr / 1e3:.4f} kN.m)")

# 3) Poteau articulé de 5 m déjà comprimé à 0,4 Pcr (cas 1) : multiplicateur du cas variable (cas 2)
L = 5.0
pcr = math.pi**2 * E * IZ / L**2
m, s, nd, b = barre(L, 1)
m.add_support(nd[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(nd[-1], dy=True, dz=True)
m.add_point_load(nd[-1], fx=-0.4 * pcr, case=1)
m.add_point_load(nd[-1], fx=-1.0, case=2)
res = m.buckling_cte_var(1, 2, n_modes=1, bisection=True)
print(f"charge variable critique : {res['factors'][0] / 1e3:.3f} kN (0,6 Pcr = {0.6 * pcr / 1e3:.3f} kN)")
