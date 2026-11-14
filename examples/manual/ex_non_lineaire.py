# Large displacements: elastic cantilever (elastica) then von Mises truss in
# snap-through, traced by arc length, and bracketing of the collapse load.
import math
import keta

# 1) Cantilever of 2 m, tip load such that P L² / (E I) = 1
L, EI = 2.0, 210e9 * 0.05**4 / 12
m = keta.Model()
nd = [m.add_node(L * i / 16, 0, 0) for i in range(17)]
sec = m.add_section_props(sx=0.05**2, iy=0.05**4 / 12, iz=0.05**4 / 12, ix=1e-6)
for i in range(16):
    m.add_beam(nd[i], nd[i + 1], sec, m.steel() if i == 0 else 1)
m.fix(nd[0])
m.add_point_load(nd[-1], fz=-EI / L**2)
nl = m.nonlinear(case=1, substeps=4)
ux, _, uz, *_ = nl.displacement(nd[-1])
print(f"elastica: uz/L = {uz / L:.5f}, ux/L = {ux / L:.5f} ({nl.iterations} iterations)")
print("reference (Bisshopp and Drucker, PL²/EI = 1): uz/L = -0.30172, ux/L = -0.05643")
lin = m.analyze().case(1).displacement(nd[-1])[2]
print(f"linear analysis: uz/L = {lin / L:.5f}")

# 2) von Mises truss: half-span a = 2 m, height h = 0.2 m, EA = 1e8 N
a, h, EA = 2.0, 0.2, 1e8
m = keta.Model()
g, s, d = m.add_node(-a, 0, 0), m.add_node(0, 0, h), m.add_node(a, 0, 0)
sec = m.add_section_props(sx=1e-3, iy=1e-8, iz=1e-8, ix=1e-8)
mat = m.add_material(e=EA / 1e-3, nu=0.3)
m.add_truss(g, s, sec, mat)
m.add_truss(s, d, sec, mat)
m.pin(g)
m.pin(d)
m.add_support(s, dy=True)
m.add_point_load(s, fz=-1000.0)                     # reference load of 1 kN
arc = m.arc_length(case=1, arc_length=0.002, max_steps=400, stop_disp=0.45)
lam = arc.lambdas()
k = next(i for i in range(1, arc.n_points - 1) if lam[i + 1] < lam[i])   # first maximum
print(arc)
print(f"limit point: lambda = {lam[k]:.3f} (point {k}), deflection {arc.node_disp(k, s)[2] * 1e3:.2f} mm")
L0 = math.hypot(a, h)
P = lambda w: -2 * EA * (math.hypot(a, h - w) - L0) / L0 * (h - w) / math.hypot(a, h - w)
Pmax = max(P(h * i / 20000) for i in range(20001))
print(f"exact equilibrium of the pinned-pinned bar: P max = {Pmax / 1000:.3f} kN")
print(f"points where the number of negative pivots changes: {arc.crossings()}; minimum of lambda: {round(min(lam), 3)}")
lo, hi, oracle = m.ruin(case=1)
print(f"collapse load bracketed: {lo:.3f} <= lambda <= {hi:.3f} ({oracle})")
