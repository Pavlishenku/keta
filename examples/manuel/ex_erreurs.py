# Détection des erreurs : identifiant absent (KeyError), donnée ou calcul invalide (ValueError).
import keta

m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(6, 0, 0)
poutre = m.add_beam(a, b, m.add_section_rect(0.2, 0.4), m.steel())
m.add_point_load(b, fz=-1e3)

try:
    m.add_section("IPE999")
except KeyError as e:
    print("KeyError   :", e)

try:
    m.add_udl(poutre, wz=-1.0, frame="oblique")
except ValueError as e:
    print("ValueError :", e)

try:
    m.analyze()                    # aucun appui : mécanisme
except ValueError as e:
    print("ValueError :", str(e).splitlines()[0])

m.fix(a)
r = m.analyze()
try:
    r.case(2)
except KeyError as e:
    print("KeyError   :", e)

m2 = keta.Model()
m2.add_node(0, 0, 0)
m2.add_beam(1, 7, 1, 1)            # références inexistantes
print("validate() :", m2.validate())
