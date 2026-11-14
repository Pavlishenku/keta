# Lignes d'influence : pont continu à trois travées 25 + 30 + 25 m, moment sur pile
# et moment à mi-travée centrale ; extrêmes sous charges ponctuelle, répartie et convois.
import keta

portees = [25.0, 30.0, 25.0]
m = keta.Model()
xs = [0.0]
for L in portees:
    xs += [xs[-1] + L * i / 10 for i in range(1, 11)]
n = [m.add_node(x, 0, 0) for x in xs]
sec = m.add_section_props(sx=6.0, iy=4.0, iz=40.0, ix=5.0)
b = [m.add_beam(n[i], n[i + 1], sec, m.concrete()) for i in range(len(n) - 1)]
m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
for i in (10, 20, 30):
    m.add_support(n[i], dy=True, dz=True)

# My à la fin de la barre 10 (sur la pile 1) et au milieu de la travée centrale (fin de la barre 15)
lig_pile = m.influence_line(b, dict(kind="bar_force", bar=b[9], x=2.5, comp=4))
lig_mi = m.influence_line(b, dict(kind="bar_force", bar=b[14], x=3.0, comp=4))
print(f"ligne : {len(lig_pile.s)} abscisses de {lig_pile.s[0]:.1f} à {lig_pile.s[-1]:.1f} m")
for nom, lig in (("M sur pile", lig_pile), ("M mi-travée", lig_mi)):
    mx, mn = lig.point_extreme()
    print(f"{nom:12}: charge unité -> max {mx:7.3f}, min {mn:7.3f} m ; "
          f"aires + {lig.integral_pos():8.3f} / - {lig.integral_neg():8.3f} m2")
    print(f"{'':12}  zones positives : " + ", ".join(f"[{a:.2f}; {c:.2f}]" for a, c, _ in lig.zones(1)))
    print(f"{'':12}  ordonnée à s = 40 m : {lig.at(40.0):.4f} m")

# Moment à mi-travée centrale : convois et charges réparties (valeurs caractéristiques, N.m)
q = 9e3 * 3.0                                    # 9 kN/m2 sur une voie de 3 m
cas = [("tandem LM1 voie 1", keta.convoy_extreme(lig_mi, keta.lm1_tandem(1), 2000)),
       ("répartie 27 kN/m", keta.udl_extreme(lig_mi, q)),
       ("LM71", keta.convoy_extreme(lig_mi, keta.rail_lm71(), 2000))]
for nom, (mx, mn, rmx, rmn) in cas:
    print(f"{nom:18}: max = {mx / 1e3:9.1f} kN.m, min = {mn / 1e3:9.1f} kN.m")
