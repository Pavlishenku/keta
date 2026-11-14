# Poutre rectangulaire 0,40 x 0,80 m en C40/50, un câble de 1500 mm² (Y1860S7) à 0,30 m sous
# le centre de gravité, tendu à 1000 MPa : moment résistant avec câble adhérent puis non
# adhérent (comparés au bloc parabole-rectangle à la main), diagramme N-M, action isostatique
# et relation moment-courbure jusqu'au point ultime.
import math
import keta

b, h, fck = 0.40, 0.80, 40e6
fcd = 0.85 * fck / 1.5
rings = [[(-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2)]]
Ap, zp, sigma_pm = 1500e-6, -0.30, 1000e6
eps_p0 = keta.tendon_initial_strain(sigma_pm)          # 1000 MPa / Ep (le béton n'est pas encore comprimé)
cable = dict(y=0.0, z=zp, Ap=Ap, eps_p0=eps_p0)

r = keta.prestressed_section_resistance(rings, fck, fcd, [cable], theta=math.pi)   # flexion positive : My < 0
print(f"adhérent : MRd = {r['M_R'] / 1e3:.1f} kN.m, sigma_p = {r['tendon_stress'][0] / 1e6:.1f} MPa, "
      f"pivot {'ABC'[min(int(r['tau']), 2)]}")
# à la main : le câble atteint la branche inclinée ; bloc 17/21 b x fcd à 99/238 x du bord comprimé
P = Ap * r["tendon_stress"][0]
x = P / (17 / 21 * b * fcd)
print(f"  bloc à la main : x = {x:.4f} m, MRd = {P * (h / 2 - zp - 99 / 238 * x) / 1e3:.1f} kN.m")

nad = dict(cable, bonded=False)                        # surtension forfaitaire de 100 MPa
r = keta.prestressed_section_resistance(rings, fck, fcd, [nad], theta=math.pi)
P = Ap * (sigma_pm + 100e6)
x = P / (17 / 21 * b * fcd)
print(f"non adhérent : MRd = {r['M_R'] / 1e3:.1f} kN.m (à la main {P * (h / 2 - zp - 99 / 238 * x) / 1e3:.1f})")

d = keta.prestressed_section_diagram(rings, fck, fcd, [cable], math.pi / 2, n_points=7)
print("diagramme (N kN, My kN.m) :", [(round(n / 1e3), round(my / 1e3)) for n, my, mz in d])

o = keta.prestress_decomposition([cable], total=(-1.5e6, 0.30e6, 0.0))
print(f"action isostatique (N, My) = ({o['isostatic'][0] / 1e3:.0f} kN, {o['isostatic'][1] / 1e3:.0f} kN.m), "
      f"part hyperstatique My = {o['hyperstatic'][1] / 1e3:.0f} kN.m")

beton = dict(kind="concrete_pr", fcd=fcd)
fibres, cables = (b, h, 1, 200), [dict(cable, A=Ap)]
ultime = keta.prestressed_ultimate_point(fibres, beton, cables, alpha=math.pi)   # fibre inférieure tendue
print(f"point ultime : kappa_u = {-ultime['kappa_y'] * 1e3:.2f} 1/km, M = {-ultime['My'] / 1e3:.1f} kN.m")
mk = keta.prestressed_m_kappa(fibres, beton, cables, [ultime["kappa_y"] * i / 5 for i in range(6)])
print("M-kappa (1/km, kN.m) :", [(round(-k * 1e3, 1), round(-m / 1e3)) for k, m, e0 in mk])
