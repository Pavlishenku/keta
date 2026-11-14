# Unilateral supports, lift-off and elastic support with a full matrix:
# two-span beam of 5 m spans, only the first span loaded.
import math
import keta

Ls, q, EI = 5.0, 10e3, 210e9 * 1e-4

def make_beam():
    m = keta.Model()
    st = m.add_material(e=210e9, nu=0.3, rho=0.0)
    s = m.add_section_props(sx=0.01, iy=1e-4, iz=1e-4, ix=1e-4)
    n = [m.add_node(Ls * i / 4, 0, 0) for i in range(9)]
    b = [m.add_beam(n[i], n[i + 1], s, st) for i in range(8)]
    m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
    m.add_support(n[4], dy=True, dz=True)
    m.add_support(n[8], dy=True, dz=True)
    for bi in b[:4]:
        m.add_udl(bi, wz=-q)
    return m, n

for option in (None, "uplift", "decol"):
    m, n = make_beam()
    if option == "uplift":
        m.set_support(n[8], uplift=dict(dz=+1))      # vertical reaction >= 0 only
    elif option == "decol":
        m.set_support(n[8], decol=("z", +1))         # the whole support is released if it lifts off
    case_ = m.analyze().case(1)
    R = dict(case_.reactions())
    Rb = R.get(n[8], (0.0,) * 6)[2]
    print(f"{str(option):6}: reactions A = {R[n[0]][2] / 1e3:7.3f}, B = {R[n[4]][2] / 1e3:7.3f}, "
          f"C = {round(Rb / 1e3, 3) + 0.0:7.3f} kN; uz(C) = {case_.displacement(n[8])[2] * 1e3:7.3f} mm")
print(f"expected free uplift: theta_B L = q L^3 / (24 EI) x L = {q * Ls**4 / (24 * EI) * 1e3:.3f} mm")

# Elastic support defined by its matrix (6 diagonal terms), rotated by 90° about Z:
# the "x" spring of the support acts along global Y
m = keta.Model()
st = m.add_material(e=210e9, nu=0.3, rho=0.0)
a, b = m.add_node(0, 0, 0), m.add_node(2, 0, 0)
m.add_truss(a, b, m.add_section_props(0.01, 1e-4, 1e-4, 1e-4), st)
m.fix(a)
m.add_support(b, dz=True, rx=True, ry=True, rz=True)
m.set_support(b, stiffness=[5e6, 0, 0, 0, 0, 0], rot=math.pi / 2)
m.add_point_load(b, fx=2e3, fy=1e3)
d = m.analyze().case(1).displacement(b)
print(f"rotated support: uy = {d[1] * 1e3:.4f} mm (1 kN / 5 MN/m = 0.2000 mm), ux = {d[0] * 1e6:.4f} µm (elongation of the bar)")
