# Fatigue: history of the midspan moment of a continuous beam 2 x 20 m under the passage
# of the FLM3 vehicle, rainflow counting and Miner damage (curve of slope m = 3).
import keta

m = keta.Model()
n = [m.add_node(2.0 * i, 0, 0) for i in range(21)]
b = [m.add_beam(n[i], n[i + 1], m.add_section_props(0.5, 0.08, 0.5, 0.05), m.steel()) for i in range(20)]
m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(n[10], dy=True, dz=True)
m.add_support(n[20], dy=True, dz=True)
line = m.influence_line(b, dict(kind="bar_force", bar=b[4], x=2.0, comp=4))  # My at x = 10 m

flm3 = keta.flm3_vehicle()
print("FLM3:", [(d, round(p / 1e3)) for d, p in flm3["axles"]], "(offset m, kN)")
r, M = keta.convoy_history(line, flm3)
print(f"history: {len(M)} positions, M from {min(M) / 1e3:.1f} to {max(M) / 1e3:.1f} kN.m")
tp = keta.turning_points(M)
print("turning points [kN.m]:", [round(x / 1e3, 1) for x in tp])
cycles = keta.rainflow(M)
for range_, mean, count in cycles:
    print(f"   cycle: range {range_ / 1e3:8.2f} kN.m, mean {mean / 1e3:8.2f} kN.m, count {count}")

W = 0.08 / 0.5                                   # section modulus I / v [m3]
dsc, NC = 71e6, 2e6                              # detail category 71 MPa at 2e6 cycles
D = sum(nb * ((e / W) / dsc) ** 3 / NC for e, _, nb in cycles)
print(f"damage of one passage: {D:.3e}; allowable passages: {1 / D:.3e}")
