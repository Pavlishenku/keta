# Fréquences propres d'une poutre sur deux appuis, nue puis portant à mi-portée un ressort,
# un appui simple ou une masse ponctuelle.
import math
import keta

E, rho, L, a = 2e11, 7800.0, 1.0, 0.05     # acier, portée, section carrée a x a
A, I, n = a * a, a**4 / 12, 80
k_ressort, masse = 1e6, 0.2                # N/m, kg

def keta_f(annexe, masse_materiau=False):
    m = keta.Model()
    acier = m.add_material(e=E, nu=0.3, rho=rho if masse_materiau else None)
    sec = m.add_section_props(sx=A, iy=I, iz=I, ix=I)
    nd = [m.add_node(L * i / n, 0, 0) for i in range(n + 1)]
    for i in range(n):
        b = m.add_beam(nd[i], nd[i + 1], sec, acier)
        if not masse_materiau:
            m.add_bar_mass(b, kind="uniform", m1=rho * A)      # translation seule : Euler-Bernoulli
    for x in nd:
        m.add_support(x, dx=True, dy=True, rx=True, rz=True)   # flexion dans le plan XZ
    m.add_support(nd[0], dz=True)
    m.add_support(nd[-1], dz=True)
    milieu = nd[n // 2]
    if annexe == "ressort":
        m.add_support(milieu, kz=k_ressort)
    elif annexe == "appui":
        m.add_support(milieu, dz=True)
    elif annexe == "masse":
        m.add_point_mass(milieu, masse)
    return m.modal(n_modes=3, method="dense").frequencies

def f_de(bl):                                # fréquence associée à beta l, l = L / 2
    return (2 * bl / L) ** 2 / (2 * math.pi) * math.sqrt(E * I / (rho * A))

def g(bl, annexe):                           # modes symétriques (demi-poutre 0 <= x <= l)
    s, c, sh, ch = math.sin(bl), math.cos(bl), math.sinh(bl), math.cosh(bl)
    if annexe == "appui":
        return s * ch - c * sh               # tan(beta l) = th(beta l)
    kr = {"ressort": k_ressort, "masse": -masse * (2 * math.pi * f_de(bl)) ** 2}.get(annexe, 0.0)
    return 4 * E * I * (2 * bl / L) ** 3 * c * ch + kr * (s * ch - c * sh)

def exactes(annexe):
    f, x, pas = [], 0.1, 1e-3
    while x < 3 * math.pi:                   # balayage en beta l, puis dichotomie sur chaque racine
        if g(x, annexe) * g(x + pas, annexe) < 0:
            u, v = x, x + pas
            for _ in range(60):
                if g(u, annexe) * g((u + v) / 2, annexe) <= 0:
                    v = (u + v) / 2
                else:
                    u = (u + v) / 2
            f.append(f_de(u))
        x += pas
    f += [f_de(j * math.pi / 2) for j in (2, 4)]      # antisymétriques : beta L = 2 pi, 4 pi
    return sorted(f)[:3]

def ligne(nom, fk, fe, ref):
    ecart = max(abs(u / v - 1) for u, v in zip(fk, fe))
    print(f"{nom:8s} keta     " + "".join(f"{v:11.4f}" for v in fk))
    print(f"{'':8s} {ref:9s}" + "".join(f"{v:11.4f}" for v in fe) + f"  écart {100 * ecart:.4f} %")
    return ecart

print(f"{'annexe':18s}" + "".join(f"{'f%d [Hz]' % j:>11s}" for j in (1, 2, 3)))
for annexe in ("aucune", "ressort", "appui", "masse"):
    assert ligne(annexe, keta_f(annexe), exactes(annexe), "exactes") < 1e-5

# Masse du matériau : l'inertie de rotation des sections est comptée (poutre de Rayleigh)
fr = [f_de(j * math.pi / 2) / math.sqrt(1 + (j * math.pi / L) ** 2 * I / A) for j in (1, 2, 3)]
assert ligne("rho", keta_f("aucune", masse_materiau=True), fr, "Rayleigh") < 5e-5
