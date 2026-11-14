# Longrine de 30 m sur sol élastique, chargée à son extrémité par une force F et un couple M.
import math
import keta

E, b, h, ks = 30e9, 0.5, 0.6, 40e6        # béton, section b x h, module de réaction du sol [N/m3]
F, M, L = 200e3, 100e3, 30.0              # force descendante, couple qui abaisse le bout, longueur
EI, k = E * b * h**3 / 12, ks * b         # k : raideur du sol par mètre de poutre [N/m2]
beta = (k / (4 * EI)) ** 0.25
print(f"beta = {beta:.5f} 1/m, beta L = {beta * L:.2f}")

m = keta.Model()
beton = m.add_material(e=E, nu=0.2)
sec = m.add_section_props(sx=b * h, iy=b * h**3 / 12, iz=h * b**3 / 12, ix=0.01)
xs = (0.0, 1.0, 2.0, 4.0, 8.0, L)
nd = [m.add_node(x, 0, 0) for x in xs]
barres = [m.add_beam(nd[i], nd[i + 1], sec, beton) for i in range(len(xs) - 1)]
for br in barres:
    m.set_soil(br, kfz=k)                 # lit de Winkler, exact quel que soit le découpage
for n in nd:
    m.add_support(n, dy=True, rx=True, rz=True)
m.add_support(nd[-1], dx=True)
m.add_point_load(nd[0], fz=-F, my=-M)
cas = m.analyze().case(1)

# Poutre semi-infinie (w compté vers le bas) : w = exp(-beta x) (C1 cos beta x + C2 sin beta x)
C1, C2 = 2 * beta * (F + beta * M) / k, -2 * beta**2 * M / k
def w(x):
    return math.exp(-beta * x) * (C1 * math.cos(beta * x) + C2 * math.sin(beta * x))
def moment(x):                            # moment, positif s'il tend la fibre supérieure
    bx = beta * x
    return 2 * beta**2 * EI * math.exp(-bx) * (C1 * math.sin(bx) - C2 * math.cos(bx))

print("   x         w keta     w exacte          My keta        My exact")
for i, x in enumerate(xs[:-1]):
    wk, Mk = -cas.displacement(nd[i])[2], cas.bar_forces(barres[i])[0][5]
    print(f"{x:4.1f} m  {wk * 1e3:8.4f} mm  {w(x) * 1e3:8.4f} mm  "
          f"{Mk / 1e3:9.3f} kN.m  {moment(x) / 1e3:9.3f} kN.m")
    assert abs(wk - w(x)) < 1e-6 * C1 and abs(Mk - moment(x)) < 1e-6 * M
ry, ry_ref = cas.displacement(nd[0])[4], -2 * beta**2 * (F + 2 * beta * M) / k
print(f"rotation de l'extrémité : keta {ry * 1e3:.5f} mrad, exacte {ry_ref * 1e3:.5f} mrad, "
      f"écart {round(100 * (ry - ry_ref) / ry_ref, 4) + 0:.4f} %")
assert abs(ry - ry_ref) < 1e-6 * abs(ry_ref)
