# Offsets: schematic composite deck, beam offset by 0.6 m below the node line.
import keta

def z(x):                     # removes the displayed negative zero
    return 0.0 if abs(x) < 1e-6 else x

L = 10.0
for e in (0.0, -0.6):
    m = keta.Model()
    a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
    beam = m.add_beam(a, b, m.add_section("HEB500"), m.steel(), ecc=(0.0, 0.0, e))
    m.add_support(a, dx=True, dy=True, dz=True, rx=True)
    m.add_support(b, dy=True, dz=True)
    m.add_point_load(b, fx=-200e3)                  # horizontal thrust on the movable support
    case_ = m.analyze().case(1)
    st = case_.bar_forces(beam)
    print(f"e = {e:5.2f} m: N = {st[5][1] / 1e3:8.2f} kN, My(mid) = {z(st[5][5]) / 1e3:8.2f} kN.m, "
          f"ux(B) = {case_.displacement(b)[0] * 1e3:7.4f} mm, ry(B) = {z(case_.displacement(b)[4]) * 1e3:7.4f} mrad")
