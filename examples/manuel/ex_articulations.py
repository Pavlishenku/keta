# Articulations d'extrémité : poutre à deux travées de 6 m, continue puis
# rendue isostatique par une rotule (ry) à l'origine de la seconde travée.
import keta

L, q = 6.0, -10e3

def z(x):                     # supprime le zéro négatif d'affichage
    return 0.0 if abs(x) < 1e-6 else x

def poutre(rotule):
    m = keta.Model()
    n = [m.add_node(L * i, 0, 0) for i in range(3)]
    sec, mat = m.add_section("IPE330"), m.steel()
    b = [m.add_beam(n[i], n[i + 1], sec, mat) for i in range(2)]
    m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
    m.add_support(n[1], dy=True, dz=True)
    m.add_support(n[2], dy=True, dz=True)
    for bi in b:
        m.add_udl(bi, wz=q)
    if rotule:
        m.set_releases(b[1], ry_or=True)   # rotation ry libérée à l'origine de la barre 2
    cas = m.analyze().case(1)
    return cas.bar_forces(b[0])[-1][5], dict(cas.reactions())[n[1]][2]

for rotule in (False, True):
    My, R = poutre(rotule)
    print(f"rotule = {rotule!s:5} : My sur appui central = {z(My) / 1e3:7.3f} kN.m, réaction = {R / 1e3:7.3f} kN")
print(f"Euler-Bernoulli : -q L²/8 = {-q * L * L / 8 / 1e3:.3f} kN.m ; réactions 10 q L/8 = {-10 * q * L / 8 / 1e3:.3f} kN et 2 q L/2 = {-q * L / 1e3:.3f} kN")
