# End hinges: two-span beam of 6 m spans, continuous, then made
# statically determinate by a hinge (ry) at the start of the second span.
import keta

L, q = 6.0, -10e3

def z(x):                     # removes the displayed negative zero
    return 0.0 if abs(x) < 1e-6 else x

def two_span(hinge):
    m = keta.Model()
    n = [m.add_node(L * i, 0, 0) for i in range(3)]
    sec, mat = m.add_section("IPE330"), m.steel()
    b = [m.add_beam(n[i], n[i + 1], sec, mat) for i in range(2)]
    m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
    m.add_support(n[1], dy=True, dz=True)
    m.add_support(n[2], dy=True, dz=True)
    for bi in b:
        m.add_udl(bi, wz=q)
    if hinge:
        m.set_releases(b[1], ry_or=True)   # rotation ry released at the start of bar 2
    case_ = m.analyze().case(1)
    return case_.bar_forces(b[0])[-1][5], dict(case_.reactions())[n[1]][2]

for hinge in (False, True):
    My, R = two_span(hinge)
    print(f"hinge = {hinge!s:5}: My over central support = {z(My) / 1e3:7.3f} kN.m, reaction = {R / 1e3:7.3f} kN")
print(f"Euler-Bernoulli: -q L²/8 = {-q * L * L / 8 / 1e3:.3f} kN.m; reactions 10 q L/8 = {-10 * q * L / 8 / 1e3:.3f} kN and 2 q L/2 = {-q * L / 1e3:.3f} kN")
