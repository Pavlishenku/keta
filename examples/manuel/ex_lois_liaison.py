# Liaisons non linéaires : suspente travaillant en traction seule, appui à jeu,
# et chemin d'équilibre d'une liaison adoucissante pilotée en déplacement.
import keta

E, I, L = 210e9, 1e-4, 3.0
kb = 3 * E * I / L**3                                  # raideur de la console en bout
for P in (-10e3, 10e3):
    m = keta.Model()
    a, c, g = m.add_node(0, 0, 0), m.add_node(L, 0, 0), m.add_node(L, 0, 1)
    m.add_beam(a, c, m.add_section_props(1e-2, I, I, 1e-4), m.add_material(E, 0.3))
    m.fix(a)
    m.fix(g)
    lien = m.add_elastic_link(g, c)                    # x du lien : de g vers c (vers le bas)
    m.set_link_law(lien, comp=0, kind="tension", k=5e6)
    m.add_point_load(c, fz=P)
    r = m.analyze_link_laws(case=1)
    print(f"P = {P / 1e3:+.0f} kN : uz = {r.displacement(c)[2] * 1e3:8.4f} mm, effort suspente = "
          f"{r.link_force(lien)[0] / 1e3:7.3f} kN ({r.iterations} itérations)")
print(f"attendus : {-10e3 / (kb + 5e6) * 1e3:.4f} mm (suspente active), {10e3 / kb * 1e3:.4f} mm (détendue)")

# Appui à jeu de 2 mm (raideur 1e8 N/m au-delà) sous la console
m = keta.Model()
a, c, g = m.add_node(0, 0, 0), m.add_node(L, 0, 0), m.add_node(L, 0, 0)
m.add_beam(a, c, m.add_section_props(1e-2, I, I, 1e-4), m.add_material(E, 0.3))
m.fix(a)
m.fix(g)
lien = m.add_elastic_link(g, c)
m.set_link_law(lien, comp=2, kind="gap", k=1e8, param=0.002)
m.add_point_load(c, fz=-20e3)
r = m.analyze_link_laws()
print(f"jeu de 2 mm : uz = {r.displacement(c)[2] * 1e3:.4f} mm, réaction du butoir {r.link_force(lien)[2] / 1e3:.3f} kN")

# Ressort adoucissant (pic à 100 N pour 10 mm, palier à 40 N) piloté en déplacement
m = keta.Model()
sol, n1 = m.add_node(0, 0, 0), m.add_node(0, 0, 0)
m.fix(sol)
m.add_support(n1, dy=True, dz=True, rx=True, ry=True, rz=True)
lien = m.add_elastic_link(sol, n1)
m.set_link_law(lien, comp=0, kind="points", points=[(0, 0), (0.01, 100.0), (0.03, 40.0), (0.1, 40.0)])
m.add_point_load(n1, fx=1.0)
chemin = m.link_laws_path(case=1, node=n1, comp=0, targets=[0.05], steps=10)
for p in chemin[::2]:
    print(f"   u = {p['u'] * 1e3:6.2f} mm : lambda = {p['lambda']:7.3f}")
