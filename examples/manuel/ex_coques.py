# Coques : dalle carrée de 6 m, simplement appuyée sur ses quatre bords, sous 10 kN/m2.
import keta

a, t, q, E, nu = 6.0, 0.20, 10e3, 30e9, 0.2
D = E * t**3 / (12 * (1 - nu**2))
print(f"solution de Navier (0,00406 q a^4 / D) : {-0.00406 * q * a**4 / D * 1e3:.4f} mm")

def dalle(n):                                      # n x n mailles, 2 triangles par maille
    m = keta.Model()
    beton = m.add_material(e=E, nu=nu)
    g = [[m.add_node(a * i / n, a * j / n, 0.0) for i in range(n + 1)] for j in range(n + 1)]
    for j in range(n):
        for i in range(n):
            n1, n2, n3, n4 = g[j][i], g[j][i + 1], g[j + 1][i + 1], g[j + 1][i]
            m.add_shell(n1, n2, n3, t, beton, slab=1)
            m.add_shell(n1, n3, n4, t, beton, slab=1)
    for j in range(n + 1):
        for i in range(n + 1):
            if i in (0, n) or j in (0, n):
                m.add_support(g[j][i], dz=True)                 # appui simple sur le contour
    m.add_support(g[0][0], dx=True, dy=True, rz=True)           # blocage des mouvements de corps rigide
    m.add_support(g[0][n], dy=True)
    m.add_surface_pressure(1, qz=-q)                            # pression verticale (repère global)

    cas = m.analyze().case(1)
    return cas.displacement(g[n // 2][n // 2])[2], cas.sum_reactions()[2]

for n in (6, 12, 24):
    w, R = dalle(n)
    print(f"{n:3d} x {n:<3d} mailles : flèche au centre = {w * 1e3:.4f} mm, réaction totale = {R / 1e3:.3f} kN")
