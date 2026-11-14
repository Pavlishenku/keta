# Flambement et second ordre : poteau bi-articulé de 6 m (HEB 200), flambement autour de l'axe faible.
# beta = pi/2 oriente l'axe faible pour que le flambement se produise selon X (voir la remarque
# sur le repère local des barres quasi verticales).
import math
import keta

L, n, P = 6.0, 12, 100e3
def poteau():
    m = keta.Model()
    nd = [m.add_node(0, 0, L * i / n) for i in range(n + 1)]
    sec, acier = m.add_section("HEB200"), m.steel()
    b = [m.add_beam(nd[i], nd[i + 1], sec, acier, beta=math.pi / 2) for i in range(n)]
    m.add_support(nd[0], dx=True, dy=True, dz=True, rz=True)
    m.add_support(nd[-1], dx=True, dy=True)
    m.add_point_load(nd[-1], fz=-P, case=1)
    return m, nd, b

m, nd, b = poteau()
fl = m.buckling(case=1, n_modes=3)
print(fl)
Iz = keta.section_props("HEB200")["Iz"]
Pcr = math.pi**2 * 210e9 * Iz / L**2
print(f"facteurs critiques : {[round(x, 4) for x in fl.factors]}")
print(f"Euler (axe faible) : Pcr / P = {Pcr / P:.4f}")
print(f"mode 1 : composantes du noeud à mi-hauteur = {[round(v, 4) + 0.0 for v in fl.mode_displacement(0, nd[n // 2])]}")

# Imperfection initiale en forme de mode 1, flèche L/300, puis calcul au second ordre
mi = m.apply_imperfection(fl, mode=0, scale=L / 300)
so = mi.second_order(case=1, max_iterations=50)
u = max(abs(so.displacement(x)[0]) for x in nd)
print(f"second ordre : {so.iterations} itérations, déplacement additionnel max {u * 1e3:.3f} mm")
amp = 1 / (1 - P / Pcr)
print(f"amplification théorique e0 / (Pcr/P - 1) = {L / 300 / (Pcr / P - 1) * 1e3:.3f} mm "
      f"(facteur 1 / (1 - P/Pcr) = {amp:.4f})")
Mz = max(abs(s[6]) for bi in b for s in so.bar_forces(bi))
print(f"moment maximal au second ordre : {Mz / 1e3:.3f} kN.m (P e0 amp = {P * L / 300 * amp / 1e3:.3f} kN.m)")
