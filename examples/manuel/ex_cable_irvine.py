# Câble pesant de 100 m (40 chaînettes, 10 kg/m, EA = 1e9 N), flèche voisine de L/100 : modes
# autour de l'équilibre (tangente de l'état en grands déplacements) contre la théorie d'Irvine.
import math
import keta

L, N, MU, EA = 100.0, 40, 10.0, 1e9
d0 = L / 100
H0 = MU * 9.80665 * L * L / (8 * d0)
l0 = L / N * (1 + 8 / 3 * (d0 / L) ** 2 - H0 / EA)        # longueur à vide d'un élément
m = keta.Model()
nd = [m.add_node(L * i / N, 0, 0) for i in range(N + 1)]
for i in range(N):
    m.add_cable(nd[i], nd[i + 1], ea=EA, w=MU * 9.80665, l0=l0)
m.fix(nd[0])
m.fix(nd[N])
for x in nd[1:N]:
    m.add_support(x, rx=True, ry=True, rz=True)
m.name_case(1, "poids propre")                           # cas vide : le poids des câbles seul

etat = m.nonlinear(case=1, substeps=20, max_iterations=80)
fleche = -etat.displacement(nd[N // 2])[2]
modes = etat.tangent().modal(n_modes=6)
H = MU * 9.80665 * L * L / (8 * fleche)                  # tension horizontale (parabole)
mu = MU * N * l0 / L                                      # masse par mètre de portée
echelle = L / math.sqrt(H / mu)
print(f"flèche {fleche:.4f} m, H = {H / 1e3:.2f} kN")
for i in range(modes.n_modes):
    dy = max(abs(d[1]) for _, d in modes.mode(i))
    dz = max(abs(d[2]) for _, d in modes.mode(i))
    print(f"mode {i + 1} : f = {modes.frequencies[i]:.4f} Hz, w L / sqrt(H/m) = {modes.omegas[i] * echelle:.3f}, "
          f"{'hors plan' if dy > dz else 'dans le plan'}")
lam2 = (mu * 9.80665 * L / H) ** 2 * L / (H * L * (1 + 8 * (fleche / L) ** 2) / EA)
lo, hi = math.pi / 2 + 1e-12, 3 * math.pi / 2 - 1e-12          # tan x = x - (4 / lambda2) x^3
for _ in range(200):
    x = 0.5 * (lo + hi)
    lo, hi = (x, hi) if math.tan(x) - x + 4 / lam2 * x ** 3 < 0 else (lo, x)
print(f"lambda2 = {lam2:.1f} ; Irvine : hors plan n pi = 3.142, 6.283, 9.425 ; antisymétrique 2 pi = 6.283 ; "
      f"symétrique {lo + hi:.3f}")
