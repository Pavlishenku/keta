# Liaisons entre noeuds : bras rigide, ressort entre deux noeuds et relation linéaire.
import keta

E = 210e9
# 1) Bras rigide : une force horizontale appliquée 0,5 m au-dessus de l'extrémité d'une console
m = keta.Model()
a, c = m.add_node(0, 0, 0), m.add_node(3, 0, 0)
m.add_beam(a, c, m.add_section_props(1e-2, 1e-4, 1e-4, 1e-4), m.add_material(E, 0.3))
m.fix(a)
s = m.add_node(3, 0, 0.5)
m.add_rigid_link(c, [s])                  # s suit le mouvement de corps rigide de c
m.add_point_load(s, fx=10e3)
uz = m.analyze().case(1).displacement(c)[2]
print(f"bras rigide : uz(c) = {uz * 1e3:.4f} mm (forme fermée {-(0.5 * 10e3) * 9 / (2 * E * 1e-4) * 1e3:.4f} mm)")

# 2) Ressort vertical entre l'extrémité d'une console et un noeud fixe
m = keta.Model()
a, c, g = m.add_node(0, 0, 0), m.add_node(3, 0, 0), m.add_node(3, 0, 0)
m.add_beam(a, c, m.add_section_props(1e-2, 1e-4, 1e-4, 1e-4), m.add_material(E, 0.3))
m.fix(a)
m.fix(g)
lien = m.add_elastic_link(g, c, kz=2e6)   # noeuds confondus : repère global
m.add_point_load(c, fz=-10e3)
cas = m.analyze().case(1)
kb = 3 * E * 1e-4 / 27
print(f"ressort : uz(c) = {cas.displacement(c)[2] * 1e3:.4f} mm (forme fermée {-10e3 / (kb + 2e6) * 1e3:.4f} mm)")
print(f"          effort dans la liaison Fz = {cas.link_force(lien)[2] / 1e3:.4f} kN")

# 3) Relation linéaire : deux consoles parallèles contraintes à la même flèche
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
cas = m.analyze().case(1)
print(f"relation : uz(c1) = {cas.displacement(c1)[2] * 1e3:.4f} mm, uz(c2) = {cas.displacement(c2)[2] * 1e3:.4f} mm")
print(f"           multiplicateur de la relation = {cas.constraint_forces()[0] / 1e3:.4f} kN")
