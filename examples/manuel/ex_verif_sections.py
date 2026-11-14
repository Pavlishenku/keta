# Outils de section : ferraillage d'une section BA en flexion composée, armatures
# de coque, moments de Wood-Armer, poutre mixte et lois uniaxiales.
import keta

# 1) Section 30 x 50 cm, deux lits d'armatures à +/- 20 cm, N = -1785 kN, My = 382 kN.m
rect = [[(-0.15, -0.25), (0.15, -0.25), (0.15, 0.25), (-0.15, 0.25)]]
lits = [(0.0, 0.20, 0.5), (0.0, -0.20, 0.5)]         # (y, z, poids) : As_i = omega x poids
d = keta.rc_design_section(rect, lits, fck=30e6, fcd=0.85 * 30e6 / 1.5, N=-1785e3, My=382e3,
                           ftd=500e6 / 1.15)
print(f"section BA : faisable = {d['feasible']}, As totale = {d['As'] * 1e4:.2f} cm2, "
      f"déformations du plan (e0, ky, kz) = {tuple(round(x, 6) + 0.0 for x in d['plane'])}")

# 2) Armatures d'un élément de coque (méthode sandwich), d = 0,16 m, bras de levier 0,15 m
a = keta.rc_shell_reinforcement(mx=40e3, my=10e3, mxy=12e3, d=0.16, z=0.15)
print("coque [cm2/m] :", {k: round(v * 1e4, 3) for k, v in a.items() if k.startswith("as")})

# 3) Moments de Wood-Armer (armatures orthogonales)
w = keta.wood_armer(mx=-22.85e3, my=24.64e3, mxy=204.53e3)
print("Wood-Armer [kN.m/m] :", {k: round(v / 1e3, 2) for k, v in w.items()})

# 4) Poutre mixte IPE 600 + dalle 1,8 m x 0,16 m, moment résistant plastique en travée
acier = dict(h=0.6, b=0.22, tw=0.012, tf=0.019, r=0.024, A=156.0e-4, Wpl=3512e-6, I=92080e-8, fy=355e6)
dalle = dict(hs=0.16, fck=30e6)
s = keta.ec4_sagging(acier, dalle, beff=1.8)
print(f"EC4 : Mpl = {s['Mpl'] / 1e3:.1f} kN.m, MRd = {s['MRd'] / 1e3:.1f} kN.m, "
      f"connexion minimale eta = {keta.ec4_eta_min(355e6, 15.0):.3f}")

# 5) Loi uniaxiale acier à écrouissage mixte sous un chargement alterné
chemin = keta.uniaxial_path("steel_mixed", dict(E=200e9, fy=355e6, H=2e9, iso=0.5),
                            [0.001, 0.003, 0.0, -0.003, 0.0])
print("acier cyclique (eps, sigma MPa) :", [(e, round(sg / 1e6, 2)) for e, sg in chemin])
