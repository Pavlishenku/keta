# Rectangular beam 0.40 x 0.80 m (C35/45), one bonded tendon of 1500 mm² and 4 HA 20 bars at the bottom
# fibre: stress limits of EN 1992-1-1 Section 7.2 and ratios, choice of the state according to f_ct,eff,
# cracking and minimum reinforcement, crack width, time-dependent losses with and without the bars.
import math

import keta

b, h, Ecm, Ep, Es = 0.40, 0.80, 34e9, 195e9, 200e9
fck, fctm, fyk, fpk = 35e6, 3.2e6, 500e6, 1860e6
rings = [[(-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2)]]
Ap, zp, As, zs = 1500e-6, -0.30, 314e-6, -0.35
tendon = dict(y=0.0, z=zp, Ap=Ap, eps_p0=1100e6 / Ep)
bars = [(y, zs, As) for y in (-0.15, -0.05, 0.05, 0.15)]
section = dict(bars=bars, Es=Es)

limits = keta.stress_limits(fck, fyk, fpk)
print(f"limits: concrete {limits['concrete_characteristic'] / 1e6:.1f} and {limits['concrete_quasi_permanent'] / 1e6:.2f} MPa, "
      f"bars {limits['bar'] / 1e6:.0f} MPa, tendons {limits['tendon'] / 1e6:.0f} MPa")

# Cracking moment (usual positive moment: theta = pi) in addition to the permanent moment Mg, without
# which the prestress alone would stress the top fibre in tension beyond f_ctm; then without prestress.
Mg = 300e3
mcr = keta.section_cracking_moment(rings, [tendon], Ecm, fctm, math.pi, sustained=(0.0, -Mg, 0.0), **section)
mrep = keta.section_cracking_moment(rings, [tendon], Ecm, fctm, math.pi, prestress=False, **section)
print(f"M_cr = {(Mg + mcr['moment']) / 1e3:.1f} kN·m, Mg included (σ_cp = {mcr['fibre']['stress'] / 1e6:.2f} MPa under P and Mg, "
      f"W = {mcr['modulus'] * 1e3:.2f} dm³); M_rep = {mrep['moment'] / 1e3:.1f} kN·m")
assert abs(mcr["moment"] - (fctm - mcr["fibre"]["stress"]) * mcr["modulus"]) < 1e-9 * mcr["moment"]

# Minimum steel against brittle failure (EN 1992-2, Expression 6.101a) with the lever arm at ULS.
z = keta.ultimate_lever_arm(rings, fck, 0.85 * fck / 1.5, [tendon], bars=bars, theta=math.pi)
print(f"lever arm of the bars at ULS {z:.3f} m, A_s,min = {keta.brittle_failure_steel(mrep['moment'], z, fyk) * 1e4:.2f} cm²")

for M in (700e3, 1000e3):
    state = keta.service_state(rings, [tendon], Ecm, fctm, sustained=(0.0, -M, 0.0), **section)
    ratios = keta.stress_ratios(state, limits)
    print(f"M = {M / 1e3:.0f} kN·m: flexural tension {state['flexural_tension'] / 1e6:.2f} MPa, "
          f"{'cracked' if state['cracked'] else 'uncracked'} state; concrete ratio {ratios['concrete_characteristic']:.2f}, "
          f"bars {ratios['bar']:.2f}, tendon {ratios['tendon']:.2f}")
    if state["cracked"]:
        w = keta.section_crack_width(rings, [tendon], Ecm, (0.0, -M, 0.0), fctm, cover=0.040, bar_dia=0.020,
                                     xi=0.5, tendon_dia=1.6 * math.sqrt(Ap), **section)
        print(f"  x = {w['x']:.3f} m, h_c,ef = {w['hc_ef'] * 1e3:.0f} mm, A_c,eff = {w['Ac_eff'] * 1e4:.0f} cm², "
              f"ρ_p,eff = {w['rho_p_eff'] * 100:.2f} %, σ_s = {w['sigma_s'] / 1e6:.0f} MPa, s_r,max = {w['s_r_max'] * 1e3:.0f} mm, "
              f"w_k = {w['wk'] * 1e3:.3f} mm")

# Time-dependent losses under the quasi-permanent moment: the bars restrain creep and shrinkage.
loss_args = dict(relaxation=[40e6], sustained=(0.0, -400e3, 0.0))
no_bars = keta.section_losses(rings, [tendon], Ecm, 2.0, -3.0e-4, **loss_args)
with_bars = keta.section_losses(rings, [tendon], Ecm, 2.0, -3.0e-4, **loss_args, **section)
s0, sz, _ = with_bars["concrete_change"]
print(f"tendon loss: {no_bars['tendons'][0]['total'] / 1e6:.1f} MPa without the bars, "
      f"{with_bars['tendons'][0]['total'] / 1e6:.1f} MPa with them; the bars take {with_bars['bars'][0]['total'] / 1e6:.1f} MPa of compression, "
      f"the concrete loses {(s0 + sz * (-h / 2)) / 1e6:.2f} MPa of compression at the bottom fibre")
