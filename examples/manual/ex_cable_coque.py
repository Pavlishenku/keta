# Tendon in a shell slab: 10 m x 2 m strip, straight tendon
# 8 cm below the mid-surface, anchored at both edges; resulting camber.
import keta

Lx, Ly, t, E = 10.0, 2.0, 0.25, 30e9
nx, ny = 20, 4
m = keta.Model()
concrete = m.add_material(e=E, nu=0.0)
g = [[m.add_node(Lx * i / nx, Ly * j / ny, 0.0) for i in range(nx + 1)] for j in range(ny + 1)]
for j in range(ny):
    for i in range(nx):
        m.add_shell(g[j][i], g[j][i + 1], g[j + 1][i + 1], t, concrete, slab=1)
        m.add_shell(g[j][i], g[j + 1][i + 1], g[j + 1][i], t, concrete, slab=1)
for j in range(ny + 1):
    m.add_support(g[j][0], dz=True)
    m.add_support(g[j][nx], dz=True)
m.add_support(g[0][0], dx=True, dy=True, rz=True)
m.add_support(g[0][nx], dy=True)

P, e = 1.0e6, -0.08
tendon = m.add_shell_tendon([(0.0, Ly / 2, e), (Lx, Ly / 2, e)], P0=P, mu=0.0, k=0.0, jack="start")
m.add_shell_prestress(tendon, factor=1.0, case=1)
case_ = m.analyze().case(1)
w = sum(case_.displacement(g[j][nx // 2])[2] for j in range(ny + 1)) / (ny + 1)
I = Ly * t**3 / 12
print(f"mean camber at midspan: {w * 1e3:.3f} mm")
print(f"equivalent beam, M = P e = {P * -e / 1e3:.0f} kN.m: M L²/(8 E I) = {P * -e * Lx**2 / (8 * E * I) * 1e3:.3f} mm")
assert abs(case_.sum_reactions()[2]) < 1e-3          # self-equilibrated equivalent loads
print("resultant of the vertical reactions: zero (self-equilibrated equivalent loads)")
