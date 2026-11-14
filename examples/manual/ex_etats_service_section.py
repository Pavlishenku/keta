# Rectangular beam 0.40 x 0.80 m (C35), one grouted bonded tendon of 1500 mm² 0.30 m below the
# centroid and 4 bars HA 20 at the bottom fibre: stresses under an increasing moment, in the
# uncracked state as long as the concrete tension stays below f_ct,eff, in the cracked state afterwards
# (EN 1992-1-1 Section 7.1(2)); the nonlinear state, concrete on the law of Section 3.1.5, for comparison.
import keta

b, h, Ecm, Ep, Es, fct_eff = 0.40, 0.80, 34e9, 195e9, 200e9, 3.2e6
rings = [[(-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2)]]
Ap, zp, As, zs = 1500e-6, -0.30, 314e-6, -0.35
tendon = dict(y=0.0, z=zp, Ap=Ap, eps_p0=1100e6 / Ep, displaced_area=Ap)  # grouted: it displaces its concrete area
bars = [(y, zs, As) for y in (-0.15, -0.05, 0.05, 0.15)]                  # (y, z, area) with omega = 1
law = dict(kind="sargin", fck=35e6)                                        # Expression (3.14), mean values

for M in (500e3, 800e3, 1100e3):                       # usual positive moments: My = -M
    section = dict(sustained=(0.0, -M, 0.0), bars=bars, Es=Es)
    state, label = keta.prestress_service_state(rings, [tendon], Ecm, **section), "uncracked"
    tension = state["fibre_max"]["stress"]
    if tension > fct_eff:
        state, label = keta.cracked_service_state(rings, [tendon], Ecm, **section), "cracked"
    print(f"M = {M / 1e3:.0f} kN·m: most tensioned fibre of the uncracked concrete {tension / 1e6:.2f} MPa, {label} state")
    print(f"  compressed concrete {state['fibre_min']['stress'] / 1e6:.2f} MPa, tendon {state['tendon_stress'][0] / 1e6:.1f} MPa, "
          f"bars {state['bar_stress'][0] / 1e6:.1f} MPa")
    nl = keta.nonlinear_service_state(rings, [tendon], law, **section)
    print(f"  law of Section 3.1.5: concrete {nl['fibre_min']['stress'] / 1e6:.2f} MPa, tendon {nl['tendon_stress'][0] / 1e6:.1f} MPa, "
          f"bars {nl['bar_stress'][0] / 1e6:.1f} MPa")
    if label == "cracked":
        # hand check: triangular block of depth x, tendon and bars, in equilibrium with M
        e0, ky, _ = state["plane"]
        x = h / 2 + e0 / ky                            # depth of the neutral axis below the top fibre
        C = 0.5 * Ecm * (e0 + ky * h / 2) * b * x      # resultant of the compressed concrete, at x / 3 from the top fibre
        Fp, Fs = Ap * state["tendon_stress"][0], 4 * As * state["bar_stress"][0]
        print(f"  neutral axis at {x:.3f} m from the top fibre; moment of the block, the tendon and the bars "
              f"{-(C * (h / 2 - x / 3) + Fp * zp + Fs * zs) / 1e3:.1f} kN·m")
        assert abs(C + Fp + Fs) < 1e-6 * abs(C) and abs(C * (h / 2 - x / 3) + Fp * zp + Fs * zs + M) < 1e-6 * M
