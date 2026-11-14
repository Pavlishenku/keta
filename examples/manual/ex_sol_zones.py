# Soil by zones and pseudo-plastic soil: quasi-rigid footing of 10 m loaded with 1 MN
# with an eccentricity, on a soil without tension, then on a soil with a capped pressure.
import keta

P, L, kz = 1e6, 10.0, 1e7                      # load [N], length [m], modulus x width [N/m per m]

def footing(e, plastic):
    m = keta.Model()
    c = m.add_material(e=2e11, nu=0.3, rho=0.0)
    s = m.add_section_props(sx=1.0, iy=1e4, iz=1e4, ix=1e4)          # quasi rigid
    xs = sorted(set([L * i / 10 for i in range(11)] + [L / 2 + e]))
    n = [m.add_node(x, 0, 0) for x in xs]
    b = [m.add_beam(n[i], n[i + 1], s, c) for i in range(len(n) - 1)]
    m.add_support(n[0], dx=True, dy=True, rx=True, rz=True)
    m.add_point_load(n[xs.index(L / 2 + e)], fz=-P, case=1)
    r = m.pse({bi: dict(zones=[dict(kfz=kz)], plastic=plastic) for bi in b}, sub_per_zone=20)
    return n, b, r.case(1)

# 1) Soil without tension (p = k w <= 0 in compression: only the 'neg' side remains), e = 3 m > L/6
n, b, case_ = footing(3.0, dict(side="neg"))
c = 3 * (L / 2 - 3.0)
print(f"no tension: {case_.iterations} iterations, yielding = {case_.plastic}")
print(f"  settlement at the loaded end {case_.displacement(n[-1])[2] * 1e3:.3f} mm "
      f"(contact over c = {c:.1f} m, maximum pressure 2P/c: {-2 * P / c / kz * 1e3:.3f} mm)")
print(f"  uplift of the opposite end {case_.displacement(n[0])[2] * 1e3:.3f} mm")
soil = case_.soil(b[-1])
print(f"  last bar: soil resultant {soil['total_global'][2] / 1e3:.2f} kN; segments (components fx, mx, fy, fz: fz):")
for s in soil["segments"][::5]:
    print(f"     [{s['x0']:.2f}; {s['x1']:.2f}] m  state {s['state'][3]:8}  mean p {s['p_mean'][3] / 1e3:9.2f} kN/m")

states = sorted({s["state"][3] for bi in b for s in case_.soil(bi)["segments"]})
print("  soil states found along the footing:", states)

# 2) Soil with a pressure capped at 160 kN/m, e = 1.5 m
n, b, case_ = footing(1.5, dict(ffz=1.6e5))
segs = [s for bi in b for s in case_.soil(bi)["segments"]]
print(f"capped: resultant {sum(case_.soil(bi)['total_global'][2] for bi in b) / 1e3:.2f} kN, "
      f"{sum(1 for s in segs if s['state'][3] == 'capped')} segments on the plateau out of {len(segs)}")
