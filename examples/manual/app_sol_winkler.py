# 30 m ground beam on an elastic soil, loaded at its end by a force F and a moment M.
import math
import keta

E, b, h, ks = 30e9, 0.5, 0.6, 40e6        # concrete, section b x h, subgrade modulus [N/m3]
F, M, L = 200e3, 100e3, 30.0              # downward force, moment pushing the end down, length
EI, k = E * b * h**3 / 12, ks * b         # k: soil stiffness per metre of beam [N/m2]
beta = (k / (4 * EI)) ** 0.25
print(f"beta = {beta:.5f} 1/m, beta L = {beta * L:.2f}")

m = keta.Model()
concrete = m.add_material(e=E, nu=0.2)
sec = m.add_section_props(sx=b * h, iy=b * h**3 / 12, iz=h * b**3 / 12, ix=0.01)
xs = (0.0, 1.0, 2.0, 4.0, 8.0, L)
nd = [m.add_node(x, 0, 0) for x in xs]
bars = [m.add_beam(nd[i], nd[i + 1], sec, concrete) for i in range(len(xs) - 1)]
for br in bars:
    m.set_soil(br, kfz=k)                 # Winkler bed, exact whatever the subdivision
for n in nd:
    m.add_support(n, dy=True, rx=True, rz=True)
m.add_support(nd[-1], dx=True)
m.add_point_load(nd[0], fz=-F, my=-M)
case = m.analyze().case(1)

# Semi-infinite beam (w positive downwards): w = exp(-beta x) (C1 cos beta x + C2 sin beta x)
C1, C2 = 2 * beta * (F + beta * M) / k, -2 * beta**2 * M / k
def w(x):
    return math.exp(-beta * x) * (C1 * math.cos(beta * x) + C2 * math.sin(beta * x))
def moment(x):                            # moment, positive when the top fibre is in tension
    bx = beta * x
    return 2 * beta**2 * EI * math.exp(-bx) * (C1 * math.sin(bx) - C2 * math.cos(bx))

print("   x         w keta      w exact          My keta        My exact")
for i, x in enumerate(xs[:-1]):
    wk, Mk = -case.displacement(nd[i])[2], case.bar_forces(bars[i])[0][5]
    print(f"{x:4.1f} m  {wk * 1e3:8.4f} mm  {w(x) * 1e3:8.4f} mm  "
          f"{Mk / 1e3:9.3f} kN.m  {moment(x) / 1e3:9.3f} kN.m")
    assert abs(wk - w(x)) < 1e-6 * C1 and abs(Mk - moment(x)) < 1e-6 * M
ry, ry_ref = case.displacement(nd[0])[4], -2 * beta**2 * (F + 2 * beta * M) / k
print(f"rotation of the end: keta {ry * 1e3:.5f} mrad, exact {ry_ref * 1e3:.5f} mrad, "
      f"diff. {round(100 * (ry - ry_ref) / ry_ref, 4) + 0:.4f} %")
assert abs(ry - ry_ref) < 1e-6 * abs(ry_ref)
