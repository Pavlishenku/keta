# Section rectangulaire en béton armé 0,30 x 0,60 m, 4 HA20 (d = 0,55 m), C30/37 :
# états I et II, moment de fissuration, ouverture des fissures, courbure, flèche,
# moment résistant, diagramme d'interaction, effort tranchant et poinçonnement.
import math
import keta

As = 4 * math.pi * 0.020**2 / 4
s = keta.RcSection(b=0.30, h=0.60, d=0.55, As=As, fck=30e6)
M = 250e3                                            # moment de service (positif en travée)
for phi in (0.0, 2.0):
    u = s.uncracked(0.0, M, phi)
    c = s.cracked(0.0, M, phi)
    print(f"phi = {phi} : Mcr = {s.cracking_moment(0.0, phi) / 1e3:7.2f} kN.m ; état I x = {u['x']:.4f} m, "
          f"I = {u['I'] * 1e4:.3f} dm4 ; état II x = {c['x']:.4f} m, sigma_c = {c['sigma_c_top'] / 1e6:6.2f} MPa, "
          f"sigma_s = {c['sigma_s'] / 1e6:6.1f} MPa")
w = s.crack_width(0.0, M, 0.0, cover=0.04, bar_dia=0.020)
print(f"ouverture des fissures : wk = {w['wk'] * 1e3:.3f} mm (sr,max = {w['s_r_max'] * 1e3:.1f} mm, "
      f"hc,ef = {w['hc_ef'] * 1e3:.1f} mm)")
k = s.curvature(0.0, M, 0.0)
print(f"courbure : kappa_I = {k['kappa_I'] * 1e3:.4f}, kappa_II = {k['kappa_II'] * 1e3:.4f}, "
      f"zeta = {k['zeta']:.4f}, kappa = {k['kappa'] * 1e3:.4f} 1/km")
L, q = 8.0, 30e3
f = s.deflection_ss(L, lambda x: q * x * (L - x) / 2, phi=2.0)
print(f"flèche à long terme d'une travée de {L} m sous {q / 1e3:.0f} kN/m : {f * 1e3:.2f} mm (L/{L / f:.0f})")
r = s.mrd(0.0)
print(f"MRd = {r['M_Rd'] / 1e3:.2f} kN.m, x = {r['x'] * 1e3:.1f} mm, aciers plastifiés : {r['steel_yields']}")
courbe = s.interaction(9)
print("interaction (N kN, MRd kN.m) :", [(round(n / 1e3), round(mm / 1e3, 1)) for n, mm in courbe])

VRd = keta.vrd_c(bw=0.30, d=0.55, Asl=As, fck=30e6)
print(f"VRd,c = {VRd / 1e3:.2f} kN ; armatures d'effort tranchant pour 300 kN (cot = 2) : "
      f"{keta.asw_s(300e3, 0.9 * 0.55, 500e6 / 1.15, 2.0) * 1e4:.3f} cm2/m ; "
      f"minimum {keta.asw_min_s(30e6, 500e6, 0.30) * 1e4:.3f} cm2/m")
nu = keta.nu_shear(30e6)
print(f"VRd,max (cot = 2) = {keta.vrd_max(0.30, 0.9 * 0.55, nu, 20e6, 2.0) / 1e3:.1f} kN")
p = keta.punching(VEd=600e3, a=0.4, b=0.4, d=0.25, rho_ly=0.01, rho_lz=0.01, fck=30e6, beta=1.15)
print(f"poinçonnement : u1 = {p['u1']:.3f} m, vEd = {p['vEd'] / 1e6:.3f} MPa, vRd,c = {p['vRd_c'] / 1e6:.3f} MPa")
