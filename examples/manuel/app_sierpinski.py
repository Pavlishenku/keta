# Treillis articulé de Sierpiński : charge au sommet et poids propre, niveaux 1 à 5.
import math
import keta

L, A, P = 24.0, 4e-3, 100e3                 # côté [m], aire des barres [m2], charge au sommet [N]

def aretes(n, a=(0.0, 0.0), b=(1.0, 0.0), c=(0.5, math.sqrt(3) / 2)):
    """Côtés des plus petits triangles du triangle de Sierpiński de niveau n (côté 1)."""
    if n == 0:
        return [(a, b), (b, c), (c, a)]
    mi = lambda p, q: ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2)
    ab, bc, ca = mi(a, b), mi(b, c), mi(c, a)
    return aretes(n - 1, a, ab, ca) + aretes(n - 1, ab, b, bc) + aretes(n - 1, ca, bc, c)

def treillis(n):
    m = keta.Model()
    acier = m.steel()
    sec = m.add_section_props(sx=A, iy=1e-8, iz=1e-8, ix=1e-8)
    nd = {}
    def noeud(p):
        cle = (round(p[0] * 2**(n + 1)), round(p[1] / (math.sqrt(3) / 2) * 2**n))   # indices de grille
        if cle not in nd:
            nd[cle] = m.add_node(p[0] * L, 0, p[1] * L)
        return nd[cle]
    barres = [m.add_truss(noeud(p), noeud(q), sec, acier) for p, q in aretes(n)]
    for k in nd.values():
        m.add_support(k, dy=True)                                 # treillis plan XZ
    gauche, droite, sommet = nd[(0, 0)], nd[(2**(n + 1), 0)], nd[(2**n, 2**n)]
    m.add_support(gauche, dx=True, dz=True)
    m.add_support(droite, dz=True)
    m.add_point_load(sommet, fz=-P, case=1)
    m.add_self_weight(case=2)
    return m, nd, barres, gauche, droite, sommet

EA = 210e9 * A
fleche = 3 * P * L / (4 * EA)               # travaux virtuels : (2/3 + 1/12) P L / EA
print(f"flèche au sommet sous P, forme fermée 3PL/4EA = {fleche * 1e3:.4f} mm")
print("niveau  noeuds  barres  actives sous P   w [mm]   poids [kN]   réactions [kN]")
for n in range(1, 6):
    m, nd, barres, gauche, droite, sommet = treillis(n)
    assert len(barres) + 3 == 2 * len(nd)                        # isostatique : b + 3 = 2 noeuds
    r = m.analyze()
    w = -r.case(1).displacement(sommet)[2]
    actives = sum(abs(r.case(1).bar_forces(b)[0][1]) > 1e-9 * P for b in barres)
    longueur = 3 * L * 1.5**n                                     # 3^(n+1) barres de longueur L / 2^n
    poids = 7850 * 9.80665 * A * longueur
    reac = dict(r.case(2).reactions())
    rg, rd = reac[gauche][2], reac[droite][2]
    print(f"{n:6d} {len(nd):7d} {len(barres):7d} {actives:15d} {w * 1e3:8.4f} {poids / 1e3:12.3f}"
          f"   {rg / 1e3:6.3f} + {rd / 1e3:6.3f}")
    assert actives == 3 * 2**n                                    # seules les barres du contour travaillent
    assert abs(w - fleche) < 1e-9 * fleche
    assert abs(rg - poids / 2) < 1e-9 * poids and abs(rd - poids / 2) < 1e-9 * poids
