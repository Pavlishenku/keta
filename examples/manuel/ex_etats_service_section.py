# Poutre rectangulaire 0,40 x 0,80 m (C35), un câble adhérent injecté de 1500 mm² à 0,30 m sous le
# centre de gravité et 4 HA 20 en fibre inférieure : contraintes sous un moment croissant, dans
# l'état non fissuré tant que la traction du béton reste sous f_ct,eff, dans l'état fissuré ensuite
# (EN 1992-1-1 § 7.1(2)) ; l'état non linéaire, béton sur la loi du § 3.1.5, à titre de comparaison.
import keta

b, h, Ecm, Ep, Es, fct_eff = 0.40, 0.80, 34e9, 195e9, 200e9, 3.2e6
rings = [[(-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2)]]
Ap, zp, As, zs = 1500e-6, -0.30, 314e-6, -0.35
cable = dict(y=0.0, z=zp, Ap=Ap, eps_p0=1100e6 / Ep, displaced_area=Ap)   # injecté : il déplace son aire de béton
barres = [(y, zs, As) for y in (-0.15, -0.05, 0.05, 0.15)]                 # (y, z, aire) avec omega = 1
loi = dict(kind="sargin", fck=35e6)                                        # expression (3.14), valeurs moyennes

for M in (500e3, 800e3, 1100e3):                       # moments positifs usuels : My = -M
    section = dict(sustained=(0.0, -M, 0.0), bars=barres, Es=Es)
    etat, nom = keta.prestress_service_state(rings, [cable], Ecm, **section), "non fissuré"
    traction = etat["fibre_max"]["stress"]
    if traction > fct_eff:
        etat, nom = keta.cracked_service_state(rings, [cable], Ecm, **section), "fissuré"
    print(f"M = {M / 1e3:.0f} kN·m : fibre la plus tendue du béton non fissuré {traction / 1e6:.2f} MPa, état {nom}")
    print(f"  béton comprimé {etat['fibre_min']['stress'] / 1e6:.2f} MPa, câble {etat['tendon_stress'][0] / 1e6:.1f} MPa, "
          f"barres {etat['bar_stress'][0] / 1e6:.1f} MPa")
    nl = keta.nonlinear_service_state(rings, [cable], loi, **section)
    print(f"  loi du § 3.1.5 : béton {nl['fibre_min']['stress'] / 1e6:.2f} MPa, câble {nl['tendon_stress'][0] / 1e6:.1f} MPa, "
          f"barres {nl['bar_stress'][0] / 1e6:.1f} MPa")
    if nom == "fissuré":
        # contrôle à la main : bloc triangulaire de hauteur x, câble et barres, en équilibre avec M
        e0, ky, _ = etat["plane"]
        x = h / 2 + e0 / ky                            # profondeur de l'axe neutre sous la fibre supérieure
        C = 0.5 * Ecm * (e0 + ky * h / 2) * b * x      # résultante du béton comprimé, à x / 3 de la fibre supérieure
        Fp, Fs = Ap * etat["tendon_stress"][0], 4 * As * etat["bar_stress"][0]
        print(f"  axe neutre à {x:.3f} m de la fibre supérieure ; moment du bloc, du câble et des barres "
              f"{-(C * (h / 2 - x / 3) + Fp * zp + Fs * zs) / 1e3:.1f} kN·m")
        assert abs(C + Fp + Fs) < 1e-6 * abs(C) and abs(C * (h / 2 - x / 3) + Fp * zp + Fs * zs + M) < 1e-6 * M
