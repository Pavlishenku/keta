# Rectangular beam 0.40 x 0.80 m in C40/50, one tendon of 1500 mm² (Y1860S7) 0.30 m below
# the centroid, stressed to 1000 MPa: resisting moment with a bonded, then an unbonded
# tendon (compared with the parabola-rectangle block by hand), N-M diagram, isostatic action
# and moment-curvature relation up to the ultimate point.
import math
import keta

b, h, fck = 0.40, 0.80, 40e6
fcd = 0.85 * fck / 1.5
rings = [[(-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2)]]
Ap, zp, sigma_pm = 1500e-6, -0.30, 1000e6
eps_p0 = keta.tendon_initial_strain(sigma_pm)          # 1000 MPa / Ep (the concrete is not yet compressed)
tendon = dict(y=0.0, z=zp, Ap=Ap, eps_p0=eps_p0)

r = keta.prestressed_section_resistance(rings, fck, fcd, [tendon], theta=math.pi)  # positive bending: My < 0
print(f"bonded: MRd = {r['M_R'] / 1e3:.1f} kN.m, sigma_p = {r['tendon_stress'][0] / 1e6:.1f} MPa, "
      f"pivot {'ABC'[min(int(r['tau']), 2)]}")
# by hand: the tendon reaches the inclined branch; block 17/21 b x fcd at 99/238 x from the compressed edge
P = Ap * r["tendon_stress"][0]
x = P / (17 / 21 * b * fcd)
print(f"  block by hand: x = {x:.4f} m, MRd = {P * (h / 2 - zp - 99 / 238 * x) / 1e3:.1f} kN.m")

unbonded = dict(tendon, bonded=False)                  # lump-sum stress increase of 100 MPa
r = keta.prestressed_section_resistance(rings, fck, fcd, [unbonded], theta=math.pi)
P = Ap * (sigma_pm + 100e6)
x = P / (17 / 21 * b * fcd)
print(f"unbonded: MRd = {r['M_R'] / 1e3:.1f} kN.m (by hand {P * (h / 2 - zp - 99 / 238 * x) / 1e3:.1f})")

d = keta.prestressed_section_diagram(rings, fck, fcd, [tendon], math.pi / 2, n_points=7)
print("diagram (N kN, My kN.m):", [(round(n / 1e3), round(my / 1e3)) for n, my, mz in d])

o = keta.prestress_decomposition([tendon], total=(-1.5e6, 0.30e6, 0.0))
print(f"isostatic action (N, My) = ({o['isostatic'][0] / 1e3:.0f} kN, {o['isostatic'][1] / 1e3:.0f} kN.m), "
      f"hyperstatic part My = {o['hyperstatic'][1] / 1e3:.0f} kN.m")

concrete = dict(kind="concrete_pr", fcd=fcd)
fibres, tendons = (b, h, 1, 200), [dict(tendon, A=Ap)]
ultimate = keta.prestressed_ultimate_point(fibres, concrete, tendons, alpha=math.pi)   # bottom fibre in tension
print(f"ultimate point: kappa_u = {-ultimate['kappa_y'] * 1e3:.2f} 1/km, M = {-ultimate['My'] / 1e3:.1f} kN.m")
mk = keta.prestressed_m_kappa(fibres, concrete, tendons, [ultimate["kappa_y"] * i / 5 for i in range(6)])
print("M-kappa (1/km, kN.m):", [(round(-k * 1e3, 1), round(-m / 1e3)) for k, m, e0 in mk])
