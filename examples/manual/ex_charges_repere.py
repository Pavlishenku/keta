# Distributed loads on an inclined bar: global, local and projected axes.
import math
import keta

alpha = math.radians(30)
L = 10.0
for frame in ("global", "local", "projected"):
    m = keta.Model()
    a = m.add_node(0, 0, 0)
    b = m.add_node(L * math.cos(alpha), 0, L * math.sin(alpha))
    bar = m.add_beam(a, b, m.add_section_rect(0.3, 0.6), m.concrete())
    m.add_support(a, dx=True, dy=True, dz=True, rx=True)
    m.add_support(b, dy=True, dz=True)
    m.add_udl(bar, wz=-10e3, frame=frame)             # 10 kN/m along -Z (or local -z)
    Fx, Fy, Fz = m.analyze().case(1).sum_reactions()
    print(f"{frame:9}: resultant of the reactions Fx = {round(Fx / 1e3, 3) + 0.0:7.3f} kN, Fz = {Fz / 1e3:7.3f} kN")
print(f"expected: global 10 L = {10 * L:.3f} kN; projected 10 L cos(alpha) = {10 * L * math.cos(alpha):.3f} kN; "
      f"local: 10 L normal to the bar")
