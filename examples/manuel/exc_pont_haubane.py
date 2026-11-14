# EXEMPLE COMPLET 2 - Pont à haubans à pylône central (structure spatiale) :
#   tablier caisson de 2 x 60 m appuyé sur les culées et sur le pylône ; pylône de 50 m
#   (20 m sous le tablier) ; 2 x 4 haubans en éventail dans le plan axial ;
#   réglage des haubans sous charges permanentes, modes propres et charge d'exploitation.
import keta

m = keta.Model()
beton, acier = m.add_material(e=36e9, nu=0.2, rho=2500.0), m.steel()
caisson = m.add_section_props(sx=6.0, iy=4.0, iz=60.0, ix=8.0, name="caisson")
fut = m.add_section_props(sx=8.0, iy=10.0, iz=10.0, ix=15.0, name="pylône")
hauban = m.add_section_props(sx=6.0e-3, iy=1e-9, iz=1e-9, ix=1e-9, name="hauban")

xs = [-60.0 + 6.0 * i for i in range(21)]
tab = {x: m.add_node(x, 0.0, 0.0) for x in xs}
barres_tablier = [m.add_beam(tab[xs[i]], tab[xs[i + 1]], caisson, beton) for i in range(20)]
pyl = [m.add_node(0.0, 0.0, z) for z in (-20.0, -10.0, 0.0, 10.0, 20.0, 30.0)]
barres_pylone = [m.add_beam(pyl[i], pyl[i + 1], fut, beton) for i in range(5)]
tete = pyl[-1]
ancrages = [-48.0, -36.0, -24.0, -12.0, 12.0, 24.0, 36.0, 48.0]
haubans = [m.add_truss(tete, tab[x], hauban, acier) for x in ancrages]

m.fix(pyl[0])                                                       # encastrement en pied de pylône
m.add_support(tab[-60.0], dy=True, dz=True, rx=True)                # culées : appuis glissants
m.add_support(tab[60.0], dy=True, dz=True, rx=True)
m.add_rigid_link(pyl[2], [tab[0.0]], dx=True, dy=True, dz=True, rx=True, ry=False, rz=False)  # appui sur le pylône

m.add_self_weight(case=1)
for b in barres_tablier:
    m.add_udl(b, wz=-60e3, case=1)                                  # superstructures
unitaires = [m.add_unit_prestress_case(h, P=1.0) for h in haubans]   # cas 2 à 9
m.name_case(1, "permanent")

# Réglage : flèche nulle sous charges permanentes aux huit points d'ancrage
cibles = [dict(kind="displacement", id=tab[x], comp=2, value=0.0) for x in ancrages]
t = m.tune_linear(base=[1], unknowns=unitaires, targets=cibles)
print("tensions de réglage (cas unitaires) [MN] :", [round(v / 1e6, 4) for v in t["x"]])
regle = m.add_combination([(1.0, 1)] + [(v, c) for v, c in zip(t["x"], unitaires)], name="permanent réglé")

# Charge d'exploitation : 3 voies (LM1 réparti 9 + 2 x 2,5 kN/m2 x 3 m) sur la travée de droite
m.add_point_load(tab[0.0], fz=0.0, case=20)
for b in barres_tablier[10:]:
    m.add_udl(b, wz=-(9e3 + 2 * 2.5e3) * 3.0, case=20)
m.name_case(20, "exploitation travée droite")
r = m.analyze(stations=5)
perm = r.combo(regle)
fmax = max(abs(perm.displacement(tab[x])[2]) for x in ancrages)
assert fmax < 1e-9
print("permanent réglé : flèches aux ancrages nulles (|uz| < 1e-9 m)")
print("efforts dans les haubans [kN] :", [round(perm.bar_forces(h)[0][1] / 1e3, 1) for h in haubans])
print(f"pylône : N en pied {perm.bar_forces(barres_pylone[0])[0][1] / 1e3:.1f} kN, "
      f"déplacement de la tête dx = {round(perm.displacement(tete)[0] * 1e3, 3) + 0.0:.3f} mm")
q = r.case(20)
print(f"exploitation : flèche maximale {min(q.displacement(n)[2] for n in tab.values()) * 1e3:.2f} mm, "
      f"dx tête de pylône {q.displacement(tete)[0] * 1e3:.2f} mm, "
      f"variation de tension max {max(abs(q.bar_forces(h)[0][1]) for h in haubans) / 1e3:.1f} kN")

mo = m.modal(n_modes=5)
for k in range(5):
    f = mo.frequencies[k]
    fr = mo.mass_ratio()[k]
    print(f"mode {k + 1} : {f:6.3f} Hz, masses effectives X {fr[0] * 100:5.1f} %  Y {fr[1] * 100:5.1f} %  "
          f"Z {fr[2] * 100:5.1f} %")
