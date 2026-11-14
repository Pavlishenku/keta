# Poutre rectangulaire 0,40 x 0,80 m (C35/45), un câble adhérent de 1500 mm² et 4 HA 20 en fibre
# inférieure : contraintes limites de l'EN 1992-1-1 § 7.2 et taux, choix de l'état selon f_ct,eff,
# fissuration et non-fragilité, ouverture de fissure, pertes différées avec et sans les barres.
import math

import keta

b, h, Ecm, Ep, Es = 0.40, 0.80, 34e9, 195e9, 200e9
fck, fctm, fyk, fpk = 35e6, 3.2e6, 500e6, 1860e6
rings = [[(-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2)]]
Ap, zp, As, zs = 1500e-6, -0.30, 314e-6, -0.35
cable = dict(y=0.0, z=zp, Ap=Ap, eps_p0=1100e6 / Ep)
barres = [(y, zs, As) for y in (-0.15, -0.05, 0.05, 0.15)]
section = dict(bars=barres, Es=Es)

limites = keta.stress_limits(fck, fyk, fpk)
print(f"limites : béton {limites['concrete_characteristic'] / 1e6:.1f} et {limites['concrete_quasi_permanent'] / 1e6:.2f} MPa, "
      f"barres {limites['bar'] / 1e6:.0f} MPa, câbles {limites['tendon'] / 1e6:.0f} MPa")

# Moment de fissuration (moment positif usuel : theta = pi) en plus du moment permanent Mg, sans
# lequel la précontrainte seule tendrait la fibre supérieure au-delà de f_ctm ; puis sans précontrainte.
Mg = 300e3
mcr = keta.section_cracking_moment(rings, [cable], Ecm, fctm, math.pi, sustained=(0.0, -Mg, 0.0), **section)
mrep = keta.section_cracking_moment(rings, [cable], Ecm, fctm, math.pi, prestress=False, **section)
print(f"M_cr = {(Mg + mcr['moment']) / 1e3:.1f} kN·m, Mg compris (σ_cp = {mcr['fibre']['stress'] / 1e6:.2f} MPa sous P et Mg, "
      f"W = {mcr['modulus'] * 1e3:.2f} dm³) ; M_rep = {mrep['moment'] / 1e3:.1f} kN·m")
assert abs(mcr["moment"] - (fctm - mcr["fibre"]["stress"]) * mcr["modulus"]) < 1e-9 * mcr["moment"]

# Acier minimal de non-fragilité (EN 1992-2, expression 6.101a) avec le bras de levier à l'ELU.
z = keta.ultimate_lever_arm(rings, fck, 0.85 * fck / 1.5, [cable], bars=barres, theta=math.pi)
print(f"bras de levier des barres à l'ELU {z:.3f} m, A_s,min = {keta.brittle_failure_steel(mrep['moment'], z, fyk) * 1e4:.2f} cm²")

for M in (700e3, 1000e3):
    etat = keta.service_state(rings, [cable], Ecm, fctm, sustained=(0.0, -M, 0.0), **section)
    taux = keta.stress_ratios(etat, limites)
    print(f"M = {M / 1e3:.0f} kN·m : traction en flexion {etat['flexural_tension'] / 1e6:.2f} MPa, "
          f"état {'fissuré' if etat['cracked'] else 'non fissuré'} ; taux béton {taux['concrete_characteristic']:.2f}, "
          f"barres {taux['bar']:.2f}, câble {taux['tendon']:.2f}")
    if etat["cracked"]:
        w = keta.section_crack_width(rings, [cable], Ecm, (0.0, -M, 0.0), fctm, cover=0.040, bar_dia=0.020,
                                     xi=0.5, tendon_dia=1.6 * math.sqrt(Ap), **section)
        print(f"  x = {w['x']:.3f} m, h_c,ef = {w['hc_ef'] * 1e3:.0f} mm, A_c,eff = {w['Ac_eff'] * 1e4:.0f} cm², "
              f"ρ_p,eff = {w['rho_p_eff'] * 100:.2f} %, σ_s = {w['sigma_s'] / 1e6:.0f} MPa, s_r,max = {w['s_r_max'] * 1e3:.0f} mm, "
              f"w_k = {w['wk'] * 1e3:.3f} mm")

# Pertes différées sous le moment quasi permanent : les barres brident le fluage et le retrait.
pertes = dict(relaxation=[40e6], sustained=(0.0, -400e3, 0.0))
seul = keta.section_losses(rings, [cable], Ecm, 2.0, -3.0e-4, **pertes)
avec = keta.section_losses(rings, [cable], Ecm, 2.0, -3.0e-4, **pertes, **section)
s0, sz, _ = avec["concrete_change"]
print(f"perte du câble : {seul['tendons'][0]['total'] / 1e6:.1f} MPa sans les barres, "
      f"{avec['tendons'][0]['total'] / 1e6:.1f} MPa avec ; les barres prennent {avec['bars'][0]['total'] / 1e6:.1f} MPa de compression, "
      f"le béton perd {(s0 + sz * (-h / 2)) / 1e6:.2f} MPa de compression en fibre inférieure")
