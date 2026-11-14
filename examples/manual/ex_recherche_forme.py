# Form finding: simply supported beam of 20 m (20 kN/m) held at midspan by a stay at 45°;
# unstressed length of the stay that cancels the midspan deflection (vertical force 5 q L / 8), then the same
# by the length defect of a prop.
import math
import keta

m = keta.Model()
mat = m.add_material(e=2e11, nu=0.3, rho=0.0)
deck = m.add_section_props(sx=0.01, iy=0.04, iz=0.04, ix=0.01)           # EA = 2e9 N
n = [m.add_node(i * 1.0, 0, 0) for i in range(21)]
for i in range(20):
    m.add_udl(m.add_beam(n[i], n[i + 1], deck, mat), wz=-2e4, case=1)
m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(n[20], dy=True, dz=True, rx=True)
pylon = m.add_node(0, 0, 10)
m.fix(pylon)
stay = m.add_cable(pylon, n[10], ea=1e8, w=0.0, h=1e5)                   # unstressed length: the unknown

r = m.tune_cables(1, [stay], [dict(id=n[10], comp=2, value=0.0)], tol=1e-10)
F = 5 * 2e4 * 20 / 8
H = F / (1 + F / 2e9)
print(f"L0 of the stay {r['x'][0]:.6f} m in {r['iterations']} iterations; midspan deflection "
      f"{r['displacement'][n[10]][2] * 1e3:.2e} mm")
print(f"by hand: vertical force 5qL/8 = {F / 1e3:.1f} kN, tension {math.hypot(H, F) / 1e3:.2f} kN")
