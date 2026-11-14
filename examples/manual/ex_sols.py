# Soil-structure interaction: earth pressure coefficients, p-y, t-z and q-z curves,
# nonlinear soil springs and earth pressure on a wall.
import math
import keta

phi = math.radians(30)
print(f"Rankine: Ka = {keta.soil_ka_rankine(phi):.4f}, Kp = {keta.soil_kp_rankine(phi):.4f}; "
      f"Jaky: K0 = {keta.soil_k0_jaky(phi):.4f}")
d = math.radians(20)
print(f"Coulomb (delta = 20°): Ka = {keta.soil_ka_coulomb(phi, d):.4f}, Kp = {keta.soil_kp_coulomb(phi, d):.4f}")
py = keta.soil_py_sand(phi=phi, gamma=10e3, X=3.0, D=1.0, k=16e6)
print("p-y sand (X = 3 m, D = 1 m), last points [(y m, p N/m)]:", [(round(y, 5), round(p)) for y, p in py[-4:]])
print("t-z clay:", [(round(z, 4), round(t)) for z, t in keta.soil_tz_clay(tmax=40e3, D=1.0)][:4], "...")
print("q-z tip:", [(round(z, 4), round(q)) for z, q in keta.soil_qz(Qp=2e6, D=1.0)][-3:])

# Pile of 10 m (tube 1.0 x 0.02) in sand: p-y springs every metre, horizontal force at the head
m = keta.Model()
nd = [m.add_node(0, 0, -1.0 * i) for i in range(11)]
tube = m.add_section_tube(d=1.0, t=0.02)
for i in range(10):
    m.add_beam(nd[i], nd[i + 1], tube, m.steel() if i == 0 else 1)
for i, x in enumerate(nd[1:], start=1):
    law = keta.soil_py_sand(phi=phi, gamma=10e3, X=float(i), D=1.0, k=16e6)
    m.add_ground_spring(x, direction=(1, 0, 0), points=[(y, p * 1.0) for y, p in law])   # 1 m of pile per spring
m.add_support(nd[-1], dz=True)
for x in nd:
    m.add_support(x, dy=True, rx=True, rz=True)
m.add_point_load(nd[0], fx=300e3)
r = m.analyze_link_laws()
print(f"pile: head displacement {r.displacement(nd[0])[0] * 1e3:.3f} mm, rotation {r.displacement(nd[0])[4] * 1e3:.4f} mrad")

# Wall of 3 m pushed 5 cm outwards: the thrust tends towards the active state
m = keta.Model()
w = [m.add_node(0, 0, -z) for z in (0.0, 1.0, 2.0, 3.0)]
for x in w:
    m.fix(x)
layer = dict(top=0.0, bottom=-3.0, sv_top=0.0, gamma=18e3, Ka=1 / 3, K0=0.5, Kp=3.0, kh_A=1e7)
links = m.add_earth_pressure(w, away=(1, 0, 0), up=(0, 0, 1), layers=[layer], distribution="lumped")
for x in w:
    m.add_support_disp(x, dx=0.05)
r = m.analyze_link_laws()
F = sum(r.link_force(l)[0] for l in links if l)
print(f"wall: total thrust {abs(F) / 1e3:.3f} kN/m (active Ka gamma H²/2 = {18e3 * 9 / 2 / 3 / 1e3:.3f} kN/m)")
