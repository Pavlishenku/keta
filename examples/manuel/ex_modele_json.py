# Commandes générales : description, validation, sauvegarde et relecture d'un modèle.
import json
import keta

m = keta.Model()
n = [m.add_node(5.0 * i, 0.0, 0.0) for i in range(3)]
s = m.add_section_rect(0.30, 0.60)
c = m.concrete()
b = [m.add_beam(n[i], n[i + 1], s, c) for i in range(2)]
for nd in n:
    m.add_support(nd, dy=True, dz=True, rx=True)
m.add_support(n[0], dx=True)
m.add_udl(b[0], wz=-20e3, case=1)
m.name_case(1, "travée 1 chargée")

print(m)                                   # résumé du modèle
print("noeuds :", m.node_count, " barres :", m.bar_count)
print("problèmes détectés :", m.validate())

texte = m.to_json()                        # modèle complet en JSON
d = json.loads(texte)
print("rubriques JSON :", sorted(k for k in d if d[k]))
m2 = keta.Model.from_json(texte)           # relecture
r1 = m.analyze().case(1).displacement(n[1])[4]
r2 = m2.analyze().case(1).displacement(n[1])[4]
print(f"rotation sur appui : {r1:.6e} rad (original), {r2:.6e} rad (relu)")
assert r1 == r2

script = m.export_script()                 # le modèle sous forme de script Python
print("script exporté :", len(script.splitlines()), "lignes ; début :")
print("\n".join(script.splitlines()[:3]))
print("version du module :", keta.__version__)
