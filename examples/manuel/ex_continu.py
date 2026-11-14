# Éléments continus : poutre-cloison en contraintes planes (quadrangles Q8)
# encastrée, soumise en bout à un effort tranchant réparti de 20 kN.
import keta

L, h, ep, P, E = 4.0, 0.4, 0.1, 20e3, 30e9
nx, ny = 20, 2
m = keta.Model()
mat = m.add_material(e=E, nu=0.0)
# noeuds au pas d'une demi-maille (les centres des mailles ne servent pas)
N = {}
for j in range(2 * ny + 1):
    for i in range(2 * nx + 1):
        if i % 2 == 0 or j % 2 == 0:
            N[i, j] = m.add_node(L * i / (2 * nx), h * j / (2 * ny), 0.0)
elements = []
for j in range(0, 2 * ny, 2):
    for i in range(0, 2 * nx, 2):
        coins = [N[i, j], N[i + 2, j], N[i + 2, j + 2], N[i, j + 2]]            # sens trigonométrique
        milieux = [N[i + 1, j], N[i + 2, j + 1], N[i + 1, j + 2], N[i, j + 1]]   # côtés 01, 12, 23, 30
        elements.append(m.add_continuum("plane_stress", coins + milieux, mat, t=ep))
for j in range(2 * ny + 1):
    m.add_support(N[0, j], dx=True, dy=True)                   # encastrement de la tranche x = 0
for e in elements[nx - 1::nx]:                                 # dernière maille de chaque rangée
    m.add_continuum_traction(e, face=1, traction=(0.0, -P / (h * ep), 0.0))   # côté x = L

cas = m.analyze().case(1)
v = cas.displacement(N[2 * nx, ny])[1]
I = ep * h**3 / 12
print(f"flèche en bout (fibre moyenne) : {v * 1e3:.4f} mm")
print(f"poutre de Timoshenko (cisaillement 5/6) : {-(P * L**3 / (3 * E * I) + P * L / (5 / 6 * E / 2 * ep * h)) * 1e3:.4f} mm")
sx, sy, txy, sz = cas.continuum_stress(elements[nx // 2], (0.0, -1.0, 0.0))  # mi-portée, bord inférieur
print(f"sigma_x à x = {L / 2 + L / (2 * nx):.1f} m, fibre inférieure : {sx / 1e6:.3f} MPa "
      f"(flexion simple : {-P * (L / 2 - L / (2 * nx)) * h / 2 / I / 1e6:.3f} MPa)")
