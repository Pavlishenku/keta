# Nonlinear links: tension-only hanger, support with a gap,
# and equilibrium path of a softening link under displacement control.
import keta

E, I, L = 210e9, 1e-4, 3.0
kb = 3 * E * I / L**3                                  # stiffness of the cantilever at its tip
for P in (-10e3, 10e3):
    m = keta.Model()
    a, c, g = m.add_node(0, 0, 0), m.add_node(L, 0, 0), m.add_node(L, 0, 1)
    m.add_beam(a, c, m.add_section_props(1e-2, I, I, 1e-4), m.add_material(E, 0.3))
    m.fix(a)
    m.fix(g)
    link = m.add_elastic_link(g, c)                    # link x: from g to c (downwards)
    m.set_link_law(link, comp=0, kind="tension", k=5e6)
    m.add_point_load(c, fz=P)
    r = m.analyze_link_laws(case=1)
    print(f"P = {P / 1e3:+.0f} kN: uz = {r.displacement(c)[2] * 1e3:8.4f} mm, hanger force = "
          f"{r.link_force(link)[0] / 1e3:7.3f} kN ({r.iterations} iterations)")
print(f"expected: {-10e3 / (kb + 5e6) * 1e3:.4f} mm (hanger active), {10e3 / kb * 1e3:.4f} mm (slack)")

# Support with a 2 mm gap (stiffness 1e8 N/m beyond it) under the cantilever
m = keta.Model()
a, c, g = m.add_node(0, 0, 0), m.add_node(L, 0, 0), m.add_node(L, 0, 0)
m.add_beam(a, c, m.add_section_props(1e-2, I, I, 1e-4), m.add_material(E, 0.3))
m.fix(a)
m.fix(g)
link = m.add_elastic_link(g, c)
m.set_link_law(link, comp=2, kind="gap", k=1e8, param=0.002)
m.add_point_load(c, fz=-20e3)
r = m.analyze_link_laws()
print(f"gap of 2 mm: uz = {r.displacement(c)[2] * 1e3:.4f} mm, reaction of the stop {r.link_force(link)[2] / 1e3:.3f} kN")

# Softening spring (peak of 100 N at 10 mm, plateau at 40 N) under displacement control
m = keta.Model()
ground, n1 = m.add_node(0, 0, 0), m.add_node(0, 0, 0)
m.fix(ground)
m.add_support(n1, dy=True, dz=True, rx=True, ry=True, rz=True)
link = m.add_elastic_link(ground, n1)
m.set_link_law(link, comp=0, kind="points", points=[(0, 0), (0.01, 100.0), (0.03, 40.0), (0.1, 40.0)])
m.add_point_load(n1, fx=1.0)
path = m.link_laws_path(case=1, node=n1, comp=0, targets=[0.05], steps=10)
for p in path[::2]:
    print(f"   u = {p['u'] * 1e3:6.2f} mm: lambda = {p['lambda']:7.3f}")
