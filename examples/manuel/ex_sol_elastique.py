# Poutre sur sol élastique (Winkler) : longrine de 40 m sous une charge centrée de 500 kN,
# comparée à la poutre infinie d'Hetényi : la solution est exacte quel que soit le maillage.
import keta

L, P, k = 40.0, 500e3, 40e6               # longueur [m], charge [N], module x largeur [N/m par m]
EI = 30e9 * 1.0 * 0.8**3 / 12
lam = (k / (4 * EI)) ** 0.25
print(f"longueur caractéristique 1/lambda = {1 / lam:.3f} m")
print(f"poutre infinie : w0 = {-P * lam / (2 * k) * 1e3:.4f} mm, M0 = {P / (4 * lam) / 1e3:.3f} kN.m")

for n in (2, 8, 40):
    m = keta.Model()
    nd = [m.add_node(L * i / n, 0, 0) for i in range(n + 1)]
    sec = m.add_section_props(sx=0.8, iy=1.0 * 0.8**3 / 12, iz=0.8 * 1.0**3 / 12, ix=0.1)  # rigide au cisaillement
    mat = m.add_material(e=30e9, nu=0.2)
    barres = [m.add_beam(nd[i], nd[i + 1], sec, mat) for i in range(n)]
    for b in barres:
        m.set_soil(b, kfz=k)              # lit élastique vertical (z local = Z)
    milieu = nd[n // 2]
    m.add_support(milieu, dx=True, dy=True, rx=True, rz=True)   # mouvements hors plan bloqués
    m.add_point_load(milieu, fz=-P)
    cas = m.analyze().case(1)
    w0 = cas.displacement(milieu)[2]
    M0 = cas.bar_forces(barres[n // 2])[0][5]
    assert abs(cas.sum_reactions()[2]) < 1e-3       # toute la charge passe dans le sol
    print(f"{n:4d} barres : w0 = {w0 * 1e3:.4f} mm, M0 = {-M0 / 1e3:.3f} kN.m")
