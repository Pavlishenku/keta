# Appuis rigides et élastiques : poutre continue de 2 x 8 m, appui central élastique.
import keta

L, q, ks = 8.0, -30e3, 20e6          # portée [m], charge [N/m], raideur [N/m]
m = keta.Model()
n = [m.add_node(L * i, 0.0, 0.0) for i in range(3)]
sec = m.add_section_props(sx=0.05, iy=2.0e-3, iz=1.0e-3, ix=1.0e-3)   # rigide à l'effort tranchant
mat = m.add_material(e=30e9, nu=0.2)
b = [m.add_beam(n[i], n[i + 1], sec, mat) for i in range(2)]

m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)   # appui fixe (rotations ry, rz libres)
m.add_support(n[1], dy=True, kz=ks)                       # appui élastique vertical
m.add_support(n[2], dy=True, dz=True)                     # appui glissant
for bi in b:
    m.add_udl(bi, wz=q)

cas = m.analyze().case(1)
for nd, (fx, fy, fz, mx, my, mz) in cas.reactions():
    print(f"noeud {nd} : Fz = {fz / 1e3:9.3f} kN")
w1 = cas.displacement(n[1])[2]
print(f"tassement de l'appui élastique : {w1 * 1e3:.3f} mm")

# Contrôle (forme fermée) : R1 = (5 q (2L)^4 / 384 EI) / (1/ks + (2L)^3 / 48 EI)
EI, Lt = 30e9 * 2.0e-3, 2 * L
R1 = (5 * -q * Lt**4 / (384 * EI)) / (1 / ks + Lt**3 / (48 * EI))
print(f"réaction attendue : {R1 / 1e3:.3f} kN")
assert abs(dict(cas.reactions())[n[1]][2] - R1) < 1e-8 * R1
assert abs(w1 + R1 / ks) < 1e-12
