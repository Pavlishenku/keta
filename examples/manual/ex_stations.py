# Result stations: subdivision of the bars, reading of the forces, stresses and
# "generalised stress" (linear combination of forces) computed in Python.
import keta

def z(x):                                         # removes the displayed negative zero
    return 0.0 if abs(x) < 1e-6 else x

L, q = 12.0, -40e3
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
beam = m.add_beam(a, b, m.add_section_rect(0.4, 1.0), m.concrete())
m.add_support(a, dx=True, dy=True, dz=True, rx=True)
m.add_support(b, dy=True, dz=True)
m.add_udl(beam, wz=q)
m.add_point_load(b, fx=-500e3)                    # compression of 500 kN

case_ = m.analyze(stations=7).case(1)             # 7 result stations per bar (step L/6)
A, I, v = 0.4 * 1.0, 0.4 * 1.0**3 / 12, 0.5
print("  x/L      x [m]      N [kN]    Vz [kN]    My [kN.m]  sig_bot [MPa]  sig_top [MPa]")
for x, N, Vy, Vz, T, My, Mz in case_.bar_forces(beam):
    s_bot = N / A + My * (-v) / I                 # sigma = N/A + My z / Iy
    s_top = N / A + My * (+v) / I
    print(f"{x:6.4f} {x * L:9.3f} {N / 1e3:11.2f} {z(Vz) / 1e3:10.2f} {z(My) / 1e3:12.2f} "
          f"{s_bot / 1e6:14.3f} {s_top / 1e6:14.3f}")
# generalised stress: C = 0.25 N + 0.33 My + 10 (any combination of forces)
C = [0.25 * s[1] + 0.33 * s[5] + 10.0 for s in case_.bar_forces(beam)]
print("generalised stress C(x):", " ".join(f"{c / 1e3:.2f}" for c in C), "(x 1000)")
print("end forces (local axes, 12 components):")
print(" ", " ".join(f"{f / 1e3:.1f}" for f in case_.end_forces(beam)))
