# Poutre rectangulaire 0,40 x 0,80 m (Ecm = 35 GPa), un câble adhérent de 1500 mm² à 0,30 m sous
# le centre de gravité, à 1200 MPa sous le poids propre : contraintes du béton aux fibres extrêmes
# et au niveau du câble, à vide puis en service, comparées à N / A + M z / I de la section
# homogénéisée écrit à la main.
import keta

b, h, Ecm, Ep = 0.40, 0.80, 35e9, 195e9
rings = [[(-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2)]]
Ac, Ic = b * h, b * h ** 3 / 12
Ap, zp, sigma_pm = 1500e-6, -0.30, 1200e6
M_vide, M_service = 250e3, 700e3                      # moments positifs usuels : My = -M (fibre inférieure tendue)

# allongement initial du câble : sa tension sous le poids propre, plus le raccourcissement du béton à son niveau
P = Ap * sigma_pm
sigma_cp = -P / Ac - P * zp * zp / Ic - M_vide * zp / Ic
cable = dict(y=0.0, z=zp, Ap=Ap, eps_p0=keta.tendon_initial_strain(sigma_pm, sigma_cp, Ecm))

# à la main : section homogénéisée (béton + n Ap au câble) sous M et la force de décompression P0 au câble
n = Ep / Ecm
Ah = Ac + n * Ap
zh = n * Ap * zp / Ah
Ih = Ic + Ac * zh ** 2 + n * Ap * (zp - zh) ** 2
P0 = Ap * Ep * cable["eps_p0"]

points = [(0.0, h / 2), (0.0, -h / 2), (0.0, zp)]     # fibre supérieure, fibre inférieure, niveau du câble
for nom, M in (("à vide", M_vide), ("en service", M_service)):
    etat = keta.prestress_service_state(rings, [cable], Ecm, sustained=(0.0, -M, 0.0), points=points)
    sup, inf, niveau = (s / 1e6 for s in etat["point_stress"])
    print(f"{nom} : fibre supérieure {sup:.2f} MPa, fibre inférieure {inf:.2f} MPa, béton au câble {niveau:.2f} MPa, "
          f"câble {etat['tendon_stress'][0] / 1e6:.1f} MPa")
    Mh = -M - P0 * zp + P0 * zh                       # moment au centre de gravité homogénéisé
    main = [-P0 / Ah + Mh * (z - zh) / Ih for _, z in points]
    print("  à la main :", ", ".join(f"{s / 1e6:.2f}" for s in main), "MPa")
    assert all(abs(s - m) < 1e-9 * abs(m) for s, m in zip(etat["point_stress"], main))
    fmax, fmin = etat["fibre_max"], etat["fibre_min"]
    print(f"  fibres extrêmes : {fmax['stress'] / 1e6:.2f} MPa à z = {fmax['z']:.2f} m, "
          f"{fmin['stress'] / 1e6:.2f} MPa à z = {fmin['z']:.2f} m")
