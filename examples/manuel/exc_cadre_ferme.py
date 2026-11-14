# EXEMPLE COMPLET 1 - Cadre fermé en béton armé sur sol élastique (ouvrage hydraulique
# ou passage inférieur), calcul plan dans XZ, tranche de 1 m.
#   ouverture entre axes 6,0 m, hauteur entre axes 4,0 m, épaisseurs 0,50 m ;
#   module de réaction du sol 30 MN/m3 ; remblai de 1,0 m sur la traverse.
import keta

B, H, e, ksol = 6.0, 4.0, 0.50, 30e6
gb, gs, K0 = 25e3, 20e3, 0.5                     # poids volumiques béton / remblai, coefficient de poussée
nb = 12                                          # barres par traverse, 8 par piédroit

m = keta.Model()
beton = m.add_material(e=32e9, nu=0.2, rho=gb / 9.80665)
sec = m.add_section_rect(b=1.0, h=e)
inf = [m.add_node(B * i / nb, 0.0, 0.0) for i in range(nb + 1)]           # radier
sup = [m.add_node(B * i / nb, 0.0, H) for i in range(nb + 1)]             # traverse
gau = [inf[0]] + [m.add_node(0.0, 0.0, H * i / 8) for i in range(1, 8)] + [sup[0]]
dro = [inf[-1]] + [m.add_node(B, 0.0, H * i / 8) for i in range(1, 8)] + [sup[-1]]
radier = [m.add_beam(inf[i], inf[i + 1], sec, beton) for i in range(nb)]
traverse = [m.add_beam(sup[i], sup[i + 1], sec, beton) for i in range(nb)]
pg = [m.add_beam(gau[i], gau[i + 1], sec, beton) for i in range(8)]
pd = [m.add_beam(dro[i], dro[i + 1], sec, beton) for i in range(8)]
for b in radier:
    m.set_soil(b, kfz=ksol * 1.0)                                          # sol sous le radier
for x in set(inf + sup + gau + dro):
    m.add_support(x, dy=True, rx=True, rz=True)                            # calcul plan XZ
m.add_support(inf[nb // 2], dx=True)                                       # blocage horizontal

# Cas 1 : poids propre ; cas 2 : remblai sur la traverse ; cas 3 : poussée des terres (au repos) ;
# cas 4 : charge d'exploitation de 10 kN/m2 sur la traverse et la poussée correspondante
m.add_self_weight(case=1)
for b in traverse:
    m.add_udl(b, wz=-gs * 1.0, case=2)
for i in range(8):
    z1, z2 = H * i / 8, H * (i + 1) / 8
    p1, p2 = K0 * gs * (1.0 + H - z1), K0 * gs * (1.0 + H - z2)            # pression à la profondeur
    m.add_linear_load(pg[i], w1=(p1, 0, 0), w2=(p2, 0, 0), case=3)
    m.add_linear_load(pd[i], w1=(-p1, 0, 0), w2=(-p2, 0, 0), case=3)
    m.add_udl(pg[i], wx=K0 * 10e3, case=4)
    m.add_udl(pd[i], wx=-K0 * 10e3, case=4)
for b in traverse:
    m.add_udl(b, wz=-10e3, case=4)
noms = ["poids propre", "remblai", "poussée des terres", "exploitation"]
for c, nom in enumerate(noms, start=1):
    m.name_case(c, nom)
elu = m.add_combination([(1.35, 1), (1.35, 2), (1.35, 3), (1.5, 4)], name="ELU")
els = m.add_combination([(1.0, 1), (1.0, 2), (1.0, 3), (1.0, 4)], name="ELS caractéristique")

r = m.analyze(stations=5)
print(m)
print("problèmes :", m.validate())
print(f"{'':22} {'M angle sup.':>13} {'M mi-trav.':>11} {'M angle inf.':>13} {'M mi-radier':>12} {'N piédroit':>11}")
kN = lambda v: round(v / 1e3, 2) + 0.0                # kN ou kN.m, sans zéro négatif
for nom, res in [(noms[c - 1], r.case(c)) for c in r.case_ids] + [("ELU", r.combo(elu)), ("ELS", r.combo(els))]:
    Ms = res.bar_forces(traverse[0])[0][5]
    Mt = res.bar_forces(traverse[nb // 2])[0][5]
    Mi = res.bar_forces(radier[0])[0][5]
    Mr = res.bar_forces(radier[nb // 2])[0][5]
    N = res.bar_forces(pg[0])[0][1]
    print(f"{nom:22} {kN(Ms):13.2f} {kN(Mt):11.2f} {kN(Mi):13.2f} {kN(Mr):12.2f} {kN(N):11.2f}")

# Pression du sol sous le radier à l'ELS : p = k w (calculée à partir des déplacements)
res = r.combo(els)
p = [-ksol * res.displacement(x)[2] / 1e3 for x in inf]
print("pression sous le radier ELS [kPa] :", " ".join(f"{v:.1f}" for v in p[::2]))
print(f"tassement moyen {sum(res.displacement(x)[2] for x in inf) / len(inf) * 1e3:.3f} mm ; "
      f"résultante des pressions {sum((p[i] + p[i + 1]) / 2 * B / nb for i in range(nb)):.2f} kN "
      f"(poids + charges : {(gb * e * (2 * B + 2 * H) + gs * B + 10e3 * B) / 1e3:.2f} kN)")
