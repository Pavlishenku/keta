# Réglage de haubans : deux haubans soutenant une poutre-console de 20 m ; on cherche
# les tensions à vide qui ramènent à zéro la flèche des deux points d'attache.
import keta

m = keta.Model()
acier = m.steel()
pied, tete = m.add_node(0, 0, 0), m.add_node(0, 0, 12)
p = [m.add_node(5.0 * i, 0, 0) for i in range(1, 5)]          # x = 5, 10, 15, 20 m
tablier = m.add_section_props(sx=0.05, iy=0.02, iz=0.2, ix=0.01)
pylone = m.add_section_props(sx=0.2, iy=0.5, iz=0.5, ix=0.3)
m.add_beam(pied, p[0], tablier, acier)
for i in range(3):
    m.add_beam(p[i], p[i + 1], tablier, acier)
m.add_beam(pied, tete, pylone, acier)
hauban = m.add_section_props(sx=4e-3, iy=1e-9, iz=1e-9, ix=1e-9)
h1 = m.add_truss(tete, p[1], hauban, acier)                    # attache à 10 m
h2 = m.add_truss(tete, p[3], hauban, acier)                    # attache à 20 m
m.fix(pied)
for b in (1, 2, 3, 4):
    m.add_udl(b, wz=-80e3, case=1)                             # poids du tablier
m.name_case(1, "permanent")
c1 = m.add_unit_prestress_case(h1, P=1.0)                      # cas unitaires de tension
c2 = m.add_unit_prestress_case(h2, P=1.0)

cibles = [dict(kind="displacement", id=p[1], comp=2, value=0.0),
          dict(kind="displacement", id=p[3], comp=2, value=0.0)]
t = m.tune_linear(base=[1], unknowns=[c1, c2], targets=cibles)
print("tensions de réglage [kN] :", [round(x / 1e3, 3) for x in t["x"]])
print("flèches obtenues [mm]    :", [round(x * 1e3, 9) + 0.0 for x in t["achieved"]])
print("flèches sans réglage [mm]:", [round(x * 1e3, 3) for x in t["base"]])

# vérification : combinaison du cas permanent et des cas unitaires pondérés
cmb = m.add_combination([(1.0, 1), (t["x"][0], c1), (t["x"][1], c2)], name="réglé")
res = m.analyze().combo(cmb)
print(f"combinaison : uz(10 m) = {res.displacement(p[1])[2] * 1e3:.6f} mm, "
      f"N hauban 1 = {res.bar_forces(h1)[0][1] / 1e3:.3f} kN, N hauban 2 = {res.bar_forces(h2)[0][1] / 1e3:.3f} kN")
