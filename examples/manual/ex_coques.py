# Shells: square slab of 6 m, simply supported on its four edges, under 10 kN/m2.
import keta

a, t, q, E, nu = 6.0, 0.20, 10e3, 30e9, 0.2
D = E * t**3 / (12 * (1 - nu**2))
print(f"Navier solution (0.00406 q a^4 / D): {-0.00406 * q * a**4 / D * 1e3:.4f} mm")

def slab(n):                                       # n x n cells, 2 triangles per cell
    m = keta.Model()
    concrete = m.add_material(e=E, nu=nu)
    g = [[m.add_node(a * i / n, a * j / n, 0.0) for i in range(n + 1)] for j in range(n + 1)]
    for j in range(n):
        for i in range(n):
            n1, n2, n3, n4 = g[j][i], g[j][i + 1], g[j + 1][i + 1], g[j + 1][i]
            m.add_shell(n1, n2, n3, t, concrete, slab=1)
            m.add_shell(n1, n3, n4, t, concrete, slab=1)
    for j in range(n + 1):
        for i in range(n + 1):
            if i in (0, n) or j in (0, n):
                m.add_support(g[j][i], dz=True)                 # simple support along the boundary
    m.add_support(g[0][0], dx=True, dy=True, rz=True)           # restraint of the rigid-body motions
    m.add_support(g[0][n], dy=True)
    m.add_surface_pressure(1, qz=-q)                            # vertical pressure (global axes)

    case_ = m.analyze().case(1)
    return case_.displacement(g[n // 2][n // 2])[2], case_.sum_reactions()[2]

for n in (6, 12, 24):
    w, R = slab(n)
    print(f"{n:3d} x {n:<3d} cells: centre deflection = {w * 1e3:.4f} mm, total reaction = {R / 1e3:.3f} kN")
