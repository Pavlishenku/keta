# Skew support: roller inclined at 30° in the XZ plane at the end of a beam.
import math
import keta

m = keta.Model()
a, c = m.add_node(0, 0, 0), m.add_node(6, 0, 0)
m.add_beam(a, c, m.add_section("IPE400"), m.steel())
m.add_support(a, dx=True, dy=True, dz=True, rx=True)
t = math.radians(30)
axes = [[math.cos(t), 0.0, math.sin(t)],        # local x: free direction of the roller
        [0.0, 1.0, 0.0],                          # local y
        [-math.sin(t), 0.0, math.cos(t)]]         # local z: normal to the rolling plane
m.add_skew_support(c, axes, dy=True, dz=True)     # restrains local y and z
m.add_point_load(c, fz=-50e3)

case_ = m.analyze().case(1)
node, loc, glob = case_.skew_reactions()[0]
print(f"node {node}: local reaction  Rz' = {loc[2] / 1e3:.3f} kN")
print(f"        global reaction Rx = {glob[0] / 1e3:.3f} kN, Rz = {glob[2] / 1e3:.3f} kN")
print(f"reaction of the fixed support: Rx = {case_.reactions()[0][1][0] / 1e3:.3f} kN")
assert abs(loc[2] - 50e3 / math.cos(t)) < 1e-6 * 50e3
