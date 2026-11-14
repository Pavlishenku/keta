# Links between nodes: rigid arm, spring between two nodes and linear constraint.
import keta

E = 210e9
# 1) Rigid arm: a horizontal force applied 0.5 m above the end of a cantilever
m = keta.Model()
a, c = m.add_node(0, 0, 0), m.add_node(3, 0, 0)
m.add_beam(a, c, m.add_section_props(1e-2, 1e-4, 1e-4, 1e-4), m.add_material(E, 0.3))
m.fix(a)
s = m.add_node(3, 0, 0.5)
m.add_rigid_link(c, [s])                  # s follows the rigid-body motion of c
m.add_point_load(s, fx=10e3)
uz = m.analyze().case(1).displacement(c)[2]
print(f"rigid arm: uz(c) = {uz * 1e3:.4f} mm (closed form {-(0.5 * 10e3) * 9 / (2 * E * 1e-4) * 1e3:.4f} mm)")

# 2) Vertical spring between the end of a cantilever and a fixed node
m = keta.Model()
a, c, g = m.add_node(0, 0, 0), m.add_node(3, 0, 0), m.add_node(3, 0, 0)
m.add_beam(a, c, m.add_section_props(1e-2, 1e-4, 1e-4, 1e-4), m.add_material(E, 0.3))
m.fix(a)
m.fix(g)
link = m.add_elastic_link(g, c, kz=2e6)   # coincident nodes: global axes
m.add_point_load(c, fz=-10e3)
case_ = m.analyze().case(1)
kb = 3 * E * 1e-4 / 27
print(f"spring: uz(c) = {case_.displacement(c)[2] * 1e3:.4f} mm (closed form {-10e3 / (kb + 2e6) * 1e3:.4f} mm)")
print(f"        force in the link Fz = {case_.link_force(link)[2] / 1e3:.4f} kN")

# 3) Linear constraint: two parallel cantilevers constrained to the same deflection
m = keta.Model()
sec, mat = m.add_section_props(1e-2, 1e-4, 1e-4, 1e-4), m.add_material(E, 0.3)
a1, c1 = m.add_node(0, 0, 0), m.add_node(3, 0, 0)
a2, c2 = m.add_node(0, 2, 0), m.add_node(3, 2, 0)
m.add_beam(a1, c1, sec, mat)
m.add_beam(a2, c2, sec, mat)
m.fix(a1)
m.fix(a2)
m.add_constraint([(c1, 2, 1.0), (c2, 2, -1.0)])      # uz(c1) - uz(c2) = 0
m.add_point_load(c1, fz=-10e3)
case_ = m.analyze().case(1)
print(f"constraint: uz(c1) = {case_.displacement(c1)[2] * 1e3:.4f} mm, uz(c2) = {case_.displacement(c2)[2] * 1e3:.4f} mm")
print(f"            multiplier of the constraint = {case_.constraint_forces()[0] / 1e3:.4f} kN")
