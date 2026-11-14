# Pin-jointed Sierpinski truss: load at the apex and self-weight, levels 1 to 5.
import math
import keta

L, A, P = 24.0, 4e-3, 100e3                 # side [m], bar area [m2], apex load [N]

def edges(n, a=(0.0, 0.0), b=(1.0, 0.0), c=(0.5, math.sqrt(3) / 2)):
    """Sides of the smallest triangles of the Sierpinski triangle of level n (side 1)."""
    if n == 0:
        return [(a, b), (b, c), (c, a)]
    mid = lambda p, q: ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2)
    ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
    return edges(n - 1, a, ab, ca) + edges(n - 1, ab, b, bc) + edges(n - 1, ca, bc, c)

def truss(n):
    m = keta.Model()
    steel = m.steel()
    sec = m.add_section_props(sx=A, iy=1e-8, iz=1e-8, ix=1e-8)
    nd = {}
    def node(p):
        key = (round(p[0] * 2**(n + 1)), round(p[1] / (math.sqrt(3) / 2) * 2**n))   # grid indices
        if key not in nd:
            nd[key] = m.add_node(p[0] * L, 0, p[1] * L)
        return nd[key]
    bars = [m.add_truss(node(p), node(q), sec, steel) for p, q in edges(n)]
    for k in nd.values():
        m.add_support(k, dy=True)                                 # plane truss XZ
    left, right, apex = nd[(0, 0)], nd[(2**(n + 1), 0)], nd[(2**n, 2**n)]
    m.add_support(left, dx=True, dz=True)
    m.add_support(right, dz=True)
    m.add_point_load(apex, fz=-P, case=1)
    m.add_self_weight(case=2)
    return m, nd, bars, left, right, apex

EA = 210e9 * A
deflection = 3 * P * L / (4 * EA)           # virtual work: (2/3 + 1/12) P L / EA
print(f"apex deflection under P, closed form 3PL/4EA = {deflection * 1e3:.4f} mm")
print("level   nodes    bars  active under P   w [mm]  weight [kN]   reactions [kN]")
for n in range(1, 6):
    m, nd, bars, left, right, apex = truss(n)
    assert len(bars) + 3 == 2 * len(nd)                          # statically determinate: b + 3 = 2 nodes
    r = m.analyze()
    w = -r.case(1).displacement(apex)[2]
    active = sum(abs(r.case(1).bar_forces(b)[0][1]) > 1e-9 * P for b in bars)
    length = 3 * L * 1.5**n                                       # 3^(n+1) bars of length L / 2^n
    weight = 7850 * 9.80665 * A * length
    reac = dict(r.case(2).reactions())
    rl, rr = reac[left][2], reac[right][2]
    print(f"{n:5d} {len(nd):7d} {len(bars):7d} {active:15d} {w * 1e3:8.4f} {weight / 1e3:12.3f}"
          f"   {rl / 1e3:6.3f} + {rr / 1e3:6.3f}")
    assert active == 3 * 2**n                                     # only the boundary bars work
    assert abs(w - deflection) < 1e-9 * deflection
    assert abs(rl - weight / 2) < 1e-9 * weight and abs(rr - weight / 2) < 1e-9 * weight
