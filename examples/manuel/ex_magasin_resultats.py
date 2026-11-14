# Portique simple calculé, enregistré, puis relu : les efforts relus sont identiques au bit près.
import os
import tempfile
import keta

m = keta.Model()
a, b, c = m.add_node(0, 0, 0), m.add_node(0, 0, 4), m.add_node(6, 0, 4)
acier = m.add_material(210e9, 0.3, rho=7850.0)
section = m.add_section_rect(0.2, 0.4)
poteau, traverse = m.add_beam(a, b, section, acier), m.add_beam(b, c, section, acier)
m.add_support(a, True, True, True, True, True, True)
m.add_support(c, False, True, True)
m.add_udl(traverse, wz=-20e3, case=1)
m.add_point_load(b, fx=15e3, case=2)

r = m.analyze(stations=5)
chemin = os.path.join(tempfile.mkdtemp(), "portique.json")
r.save(chemin, note="portique, variante A")
info = keta.results_info(chemin)
print(f"schéma {info['schema_version']}, note « {info['note']} », {info['cases']} cas")

relu = m.read_results(chemin)
identique = all(relu.case(i).bar_forces(barre) == r.case(i).bar_forces(barre)
                for i in (1, 2) for barre in (poteau, traverse))
print("efforts relus identiques :", identique)
print(f"moment en tête de poteau, cas 1 : {relu.case(1).bar_forces(poteau)[-1][5] / 1e3:.2f} kN.m")
