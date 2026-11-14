# Câble dans une dalle en coques : bande de 10 m x 2 m, câble rectiligne
# 8 cm sous le feuillet moyen, ancré aux deux bords ; contre-flèche obtenue.
import keta

Lx, Ly, t, E = 10.0, 2.0, 0.25, 30e9
nx, ny = 20, 4
m = keta.Model()
beton = m.add_material(e=E, nu=0.0)
g = [[m.add_node(Lx * i / nx, Ly * j / ny, 0.0) for i in range(nx + 1)] for j in range(ny + 1)]
for j in range(ny):
    for i in range(nx):
        m.add_shell(g[j][i], g[j][i + 1], g[j + 1][i + 1], t, beton, slab=1)
        m.add_shell(g[j][i], g[j + 1][i + 1], g[j + 1][i], t, beton, slab=1)
for j in range(ny + 1):
    m.add_support(g[j][0], dz=True)
    m.add_support(g[j][nx], dz=True)
m.add_support(g[0][0], dx=True, dy=True, rz=True)
m.add_support(g[0][nx], dy=True)

P, e = 1.0e6, -0.08
cable = m.add_shell_tendon([(0.0, Ly / 2, e), (Lx, Ly / 2, e)], P0=P, mu=0.0, k=0.0, jack="start")
m.add_shell_prestress(cable, factor=1.0, case=1)
cas = m.analyze().case(1)
w = sum(cas.displacement(g[j][nx // 2])[2] for j in range(ny + 1)) / (ny + 1)
I = Ly * t**3 / 12
print(f"contre-flèche moyenne à mi-portée : {w * 1e3:.3f} mm")
print(f"poutre équivalente, M = P e = {P * -e / 1e3:.0f} kN.m : M L²/(8 E I) = {P * -e * Lx**2 / (8 * E * I) * 1e3:.3f} mm")
assert abs(cas.sum_reactions()[2]) < 1e-3          # charges équivalentes autoéquilibrées
print("résultante des réactions verticales : nulle (charges équivalentes autoéquilibrées)")
