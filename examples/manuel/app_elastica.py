# Élastica : mât encastré en pied, libre en tête, comprimé au-delà de sa charge critique.
import math
import keta

E, a, L, n = 210e9, 0.02, 2.5, 80          # acier, section carrée a x a, hauteur, nombre de barres
EI = E * a**4 / 12
Fcr = math.pi**2 * EI / (4 * L**2)

def integrales(p):
    """Intégrales elliptiques complètes K(p) et E(p) (moyenne arithmético-géométrique)."""
    x, y, s, puissance = 1.0, math.sqrt(1 - p * p), p * p / 2, 1.0
    while x - y > 1e-15:
        c = (x - y) / 2
        x, y = (x + y) / 2, math.sqrt(x * y)
        s += puissance * c * c
        puissance *= 2
    K = math.pi / (2 * x)
    return K, K * (1 - s)

m = keta.Model()
sec = m.add_section_props(sx=a * a, iy=a**4 / 12, iz=a**4 / 12, ix=0.141 * a**4)
acier = m.add_material(e=E, nu=0.3)
nd = [m.add_node(0, 0, L * i / n) for i in range(n + 1)]
for i in range(n):
    m.add_beam(nd[i], nd[i + 1], sec, acier)
m.fix(nd[0])
for x in nd[1:]:
    m.add_support(x, dy=True, rx=True, rz=True)              # flexion dans le plan XZ
forme = [(1 - math.cos(math.pi * i / (2 * n)), 0, 0, 0, 0, 0) for i in range(n + 1)]   # mode 1
m.add_point_load(nd[-1], fz=-1.0, case=2)            # cas 2 : charge unité, pour le flambement
Fk = m.buckling(case=2, n_modes=1).factors[0]
print(f"charge critique : keta {Fk:.2f} N, pi2 EI / 4 L2 = {Fcr:.2f} N")
print("angle  F / Fcr   angle keta    dX keta   dX exact    dZ keta   dZ exact")
ecarts = []
for alpha in (30, 60, 90, 120):                              # angle visé de la tête [degrés]
    p = math.sin(math.radians(alpha) / 2)
    K, Ep = integrales(p)
    k = K / L                                                # k L = K(p), avec k = sqrt(F / EI)
    F = EI * k * k
    mi = m.apply_imperfection_shape(forme, nd, scale=L * 1e-5)   # amorce : L / 100 000 en tête
    mi.add_point_load(nd[-1], fz=-F, case=1)
    ux, _, uz, _, ry, _ = mi.nonlinear(case=1, substeps=50, max_iterations=200).displacement(nd[-1])
    dX, dZ = 2 * p / k, 2 * Ep / k - 2 * L
    print(f"{alpha:4d} °  {F / Fcr:7.4f}  {math.degrees(ry):8.3f} °  {ux:8.4f} m  {dX:8.4f} m  "
          f"{uz:8.4f} m  {dZ:8.4f} m")
    ecarts += [abs(math.degrees(ry) / alpha - 1), abs(ux / dX - 1), abs(uz / dZ - 1)]
print(f"écart relatif maximal : {100 * max(ecarts):.3f} %")
assert max(ecarts) < 1e-3
