# Poutre rectangulaire 0,40 x 0,80 m en C40/50 (Ecm = 35 GPa), deux lits de torons prétendus à
# 1400 MPa sur le banc : état au relâchement (raccourcissement élastique), pertes différées à
# 90, 365 et 20 000 jours, puis moment résistant avant et après pertes.
import math
import keta

b, h, fck, Ecm, Ep = 0.40, 0.80, 40e6, 35e9, 195e9
fcd = 0.85 * fck / 1.5
rings = [[(-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2)]]
banc = 1400e6                                             # tension sur le banc : eps_p0 = banc / Ep
lits = [dict(y=0.0, z=-0.32, Ap=900e-6, eps_p0=banc / Ep), dict(y=0.0, z=-0.25, Ap=600e-6, eps_p0=banc / Ep)]

etat = keta.prestress_service_state(rings, lits, Ecm)
for lit, sp, sc in zip(lits, etat["tendon_stress"], etat["concrete_stress"]):
    print(f"z = {lit['z']:.2f} m : torons {sp / 1e6:.1f} MPa (perte élastique {(banc - sp) / 1e6:.1f}), "
          f"béton {sc / 1e6:.2f} MPa")
# à la main pour un seul lit d'aire totale au barycentre des torons : sigma = banc / (1 + n Ap (1/Ac + e²/Ic))
Ap, e = 1500e-6, (900 * 0.32 + 600 * 0.25) / 1500
print(f"  lit unique équivalent : {banc / (1 + Ep / Ecm * Ap * (1 / (b * h) + e * e / (b * h ** 3 / 12))) / 1e6:.1f} MPa")

histoire = keta.section_loss_history(rings, lits, [90.0, 365.0, 20000.0], fck, rh=70.0, h0=0.267, t0=28.0, Ecm=Ecm)
for age, pertes in histoire:
    print(f"{age:7.0f} j :", ", ".join(
        f"{p['total'] / 1e6:6.1f} MPa (fluage {p['creep'] / 1e6:.1f}, retrait {p['shrinkage'] / 1e6:.1f}, "
        f"relaxation {p['relaxation'] / 1e6:.1f})" for p in pertes))

finales = histoire[-1][1]
apres = [dict(lit, eps_p0=p["eps_p0_after"]) for lit, p in zip(lits, finales)]
for nom, cables in (("avant pertes", lits), ("après pertes", apres)):
    r = keta.prestressed_section_resistance(rings, fck, fcd, cables, theta=math.pi)
    print(f"MRd {nom} : {r['M_R'] / 1e3:.1f} kN.m")
print(f"post-tension, 4 câbles tendus l'un après l'autre : part moyenne {keta.sequential_tensioning_factor(4):.3f}")
