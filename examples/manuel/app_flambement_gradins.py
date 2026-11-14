# Flambement d'une colonne télescopique encastrée en pied, libre en tête, de 1 à 8 tronçons de
# section décroissante, comprimée par une force verticale en tête.
import math
import keta

E = 206e9
TRONCONS = [  # (longueur [m], aire [m2], inertie 1 [m4], inertie 2 [m4]), du pied vers la tête
    (15.3, 0.0906, 0.04326, 0.03466), (15.5, 0.0645, 0.02609, 0.02281),
    (15.4, 0.0529, 0.01906, 0.01703), (15.6, 0.0478, 0.01441, 0.01335),
    (15.4, 0.0415, 0.01039, 0.01043), (15.3, 0.0346, 0.00749, 0.00756),
    (15.1, 0.0296, 0.00541, 0.00560), (14.8, 0.0302, 0.00463, 0.00489)]

def keta_pcr(n, decoupage=4):
    m = keta.Model()
    acier = m.add_material(e=E, nu=0.3)
    pied = m.add_node(0, 0, 0)
    m.fix(pied)
    z, prec = 0.0, pied
    for L, A, I1, I2 in TRONCONS[:n]:
        sec = m.add_section_props(sx=A, iy=I1, iz=I2, ix=I1 + I2)
        for _ in range(decoupage):                  # 4 barres par tronçon
            z += L / decoupage
            nn = m.add_node(0, 0, z)
            m.add_beam(prec, nn, sec, acier)
            prec = nn
    m.add_point_load(prec, fz=-1e3)                 # 1 kN : le facteur critique est Pcr en kN
    return m.buckling(n_modes=1).factors[0]

def residu(P, n, plan):
    # u = déport de la tête - flèche : u'' + (P / EI) u = 0 dans chaque tronçon,
    # u = 1 et u' = 0 au pied (encastrement), u = 0 en tête (point d'application de P)
    u, du = 1.0, 0.0
    for t in TRONCONS[:n]:
        L, k = t[0], math.sqrt(P / (E * t[2 + plan]))
        c, s = math.cos(k * L), math.sin(k * L)
        u, du = u * c + du / k * s, -u * k * s + du * c
    return u

def pcr_exacte(n):
    valeurs = []
    for plan in (0, 1):                             # chacun des deux plans principaux
        P = 1e3
        while residu(P * 1.01, n, plan) > 0:        # première racine : premier changement de signe
            P *= 1.01
        a, b = P, P * 1.01                          # puis dichotomie
        for _ in range(60):
            c = (a + b) / 2
            if residu(c, n, plan) > 0:
                a = c
            else:
                b = c
        valeurs.append(a)
    return min(valeurs) / 1e3

L1, I2 = TRONCONS[0][0], TRONCONS[0][3]                # inertie la plus faible du tronçon 1
print(f"1 tronçon, Euler pi2 E I / 4 L2 = {math.pi**2 * E * I2 / (4 * L1**2) / 1e3:.1f} kN")
print(" n   hauteur     Pcr keta      Pcr exacte     écart")
for n in range(1, 9):
    pk, pe = keta_pcr(n), pcr_exacte(n)
    ecart = 100 * (pk - pe) / pe
    H = sum(t[0] for t in TRONCONS[:n])
    print(f"{n:2d}   {H:6.1f} m  {pk:10.1f} kN  {pe:10.1f} kN  {round(ecart, 4) + 0:7.4f} %")
    assert abs(ecart) < 0.01
