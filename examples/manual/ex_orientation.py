# Section orientation (angle beta): the same IPE 300 cantilever bent
# about its strong axis (beta = 0), then about its weak axis (beta = pi/2).
import math
import keta

def z(x):                     # removes the displayed negative zero
    return 0.0 if abs(x) < 1e-6 else x

for beta in (0.0, math.pi / 2):
    m = keta.Model()
    a, b = m.add_node(0, 0, 0), m.add_node(3, 0, 0)
    m.add_beam(a, b, m.add_section("IPE300"), m.steel(), beta=beta)
    m.fix(a)
    m.add_point_load(b, fz=-5e3)
    case_ = m.analyze().case(1)
    uz = case_.displacement(b)[2]
    x, N, Vy, Vz, T, My, Mz = case_.bar_forces(1)[0]
    print(f"beta = {beta:.4f} rad: uz = {uz * 1e3:8.3f} mm, My(0) = {z(My) / 1e3:7.3f} kN.m, Mz(0) = {z(Mz) / 1e3:7.3f} kN.m")

# Local axes of a vertical bar and of an inclined bar (read through the deflection under a load)
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(0, 0, 4)
m.add_beam(a, b, m.add_section_rect(0.3, 0.5), m.concrete())
m.fix(a)
m.add_point_load(b, fx=10e3)          # horizontal force at the column head
x, N, Vy, Vz, T, My, Mz = m.analyze().case(1).bar_forces(1)[0]
print(f"column, base: Vy = {z(Vy) / 1e3:.3f} kN, Vz = {Vz / 1e3:.3f} kN, My = {My / 1e3:.3f} kN.m, Mz = {z(Mz) / 1e3:.3f} kN.m")
