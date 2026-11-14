# Rigid and elastic supports: continuous beam of 2 x 8 m, elastic central support.
import keta

L, q, ks = 8.0, -30e3, 20e6          # span [m], load [N/m], stiffness [N/m]
m = keta.Model()
n = [m.add_node(L * i, 0.0, 0.0) for i in range(3)]
sec = m.add_section_props(sx=0.05, iy=2.0e-3, iz=1.0e-3, ix=1.0e-3)   # rigid in shear
mat = m.add_material(e=30e9, nu=0.2)
b = [m.add_beam(n[i], n[i + 1], sec, mat) for i in range(2)]

m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)   # pinned support (rotations ry, rz free)
m.add_support(n[1], dy=True, kz=ks)                       # elastic vertical support
m.add_support(n[2], dy=True, dz=True)                     # sliding support
for bi in b:
    m.add_udl(bi, wz=q)

case_ = m.analyze().case(1)
for nd, (fx, fy, fz, mx, my, mz) in case_.reactions():
    print(f"node {nd}: Fz = {fz / 1e3:9.3f} kN")
w1 = case_.displacement(n[1])[2]
print(f"settlement of the elastic support: {w1 * 1e3:.3f} mm")

# Check (closed form): R1 = (5 q (2L)^4 / 384 EI) / (1/ks + (2L)^3 / 48 EI)
EI, Lt = 30e9 * 2.0e-3, 2 * L
R1 = (5 * -q * Lt**4 / (384 * EI)) / (1 / ks + Lt**3 / (48 * EI))
print(f"expected reaction: {R1 / 1e3:.3f} kN")
assert abs(dict(case_.reactions())[n[1]][2] - R1) < 1e-8 * R1
assert abs(w1 + R1 / ks) < 1e-12
