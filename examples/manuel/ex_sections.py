# Sections : catalogue, sections paramétriques et caractéristiques directes.
import json
import keta

for nom in ("IPE300", "HEA400", "HEB1000"):
    p = keta.section_props(nom)
    print(f"{nom:8} A = {p['A'] * 1e4:8.2f} cm2  Iy = {p['Iy'] * 1e8:10.1f} cm4  "
          f"Iz = {p['Iz'] * 1e8:9.1f} cm4  J = {p['J'] * 1e8:7.2f} cm4")

m = keta.Model()
acier = m.steel()
ids = {
    "catalogue HEB 300": m.add_section("heb 300"),          # casse et espaces tolérés
    "rectangle 0,4 x 0,8": m.add_section_rect(b=0.4, h=0.8),
    "tube 323,9 x 10": m.add_section_tube(d=0.3239, t=0.010),
    "I soudé 1200 x 400": m.add_section_i(h=1.2, b=0.4, tw=0.012, tf=0.030),
    "caractéristiques directes": m.add_section_props(sx=0.25, iy=0.02, iz=0.01, ix=0.015, sy=0.2, sz=0.2,
                                                    name="caisson"),
}
# les aires de cisaillement des sections paramétriques sont fixées à la création
# de la première barre (elles dépendent du coefficient de Poisson du matériau)
n0 = m.add_node(0, 0, 0)
for i, s in enumerate(ids.values()):
    m.add_beam(n0, m.add_node(1.0, i, 0), s, acier)
props = {s["id"]: s["props"] for s in json.loads(m.to_json())["sections"]}
for nom, s in ids.items():
    p = props[s]
    print(f"{nom:27} A = {p['A']:.5e}  Iy = {p['Iy']:.5e}  Avz = {p['Avz']:.5e}  J = {p['J']:.5e}")
print("convention des aires de cisaillement :", keta.shear_areas())
