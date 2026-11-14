# Road live loads on an influence line: LM1 load model per lane,
# A(l), Bc and Mc120 systems and Courbon's transverse distribution.
import keta

# Simply supported span of 30 m: midspan moment (exact triangular line, sampled)
L = 30.0
m = keta.Model()
n = [m.add_node(L * i / 12, 0, 0) for i in range(13)]
b = [m.add_beam(n[i], n[i + 1], m.add_section_props(5.0, 3.0, 30.0, 4.0), m.concrete()) for i in range(12)]
m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(n[-1], dy=True, dz=True)
line = m.influence_line(b, dict(kind="bar_force", bar=b[5], x=2.5, comp=4), per_bar=21)
print(f"maximum ordinate |eta| = {abs(line.point_extreme()[1]):.4f} m (L/4 = {L / 4:.4f} m)")

# EN 1991-2, LM1: lanes of 3 m on a carriageway of 9 m
nl, width, rest = keta.en_lane_division(9.0)
print(f"carriageway of 9 m: {nl} lanes of {width} m, remaining area {rest} m")
for lane in (1, 2, 3):
    res = keta.lm1_lane_extreme(line, lane, width)
    print(f"lane {lane}: tandem {keta.lm1_tandem(lane)['axles'][0][1] / 1e3:.0f} kN/axle, "
          f"UDL {keta.lm1_udl(lane, width) / 1e3:.1f} kN/m -> M = {min(res[:2]) / 1e3:9.1f} kN.m")

# Fascicule 61 Titre II: A(l) on the most unfavourable loaded length, Bc, Mc120
density = lambda Lc: keta.fasc61_a_lane_load(Lc, 1, 3.5) * 2   # 2 lanes, first-class bridge [N/m]
val, Lc, zones = keta.intensity_udl_extreme(line, density, sign=-1)
print(f"A(l): loaded length {Lc:.2f} m, A(30) = {keta.fasc61_a_of_l(30.0):.1f} N/m2 -> M = {val / 1e3:.1f} kN.m")
G, S = 30 * 150e3, 2 * 60 * 9.80665e3                            # self-weight and maximum Bc load
delta = keta.fasc61_delta(L, G, S)
mx, mn, rmx, rmn = keta.convoy_extreme(line, keta.fasc61_bc(), 3000)
print(f"Bc (one line of trucks): M = {mn / 1e3:.1f} kN.m, head of the convoy at s = {rmn:.2f} m; "
      f"dynamic amplification delta = {delta:.4f}")
mx, mn, rmx, rmn = keta.convoy_extreme(line, keta.fasc61_mc120(), 3000)
print(f"Mc120: M = {mn / 1e3:.1f} kN.m")

# Courbon's transverse distribution: 4 girders, unit load at y = 3 m
K = keta.courbon_coefficients([-4.5, -1.5, 1.5, 4.5], y=3.0)
print("Courbon coefficients:", [round(k, 4) for k in K])
