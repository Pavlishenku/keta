# Rectangular beam 0.40 x 0.80 m (Ecm = 35 GPa), one bonded tendon of 1500 mm² 0.30 m below
# the centroid, at 1200 MPa under self-weight: concrete stresses at the extreme fibres
# and at the tendon level, unloaded then in service, compared with N / A + M z / I of the
# transformed section written by hand.
import keta

b, h, Ecm, Ep = 0.40, 0.80, 35e9, 195e9
rings = [[(-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2)]]
Ac, Ic = b * h, b * h ** 3 / 12
Ap, zp, sigma_pm = 1500e-6, -0.30, 1200e6
M_unloaded, M_service = 250e3, 700e3                  # usual positive moments: My = -M (bottom fibre in tension)

# initial strain of the tendon: its tension under self-weight, plus the concrete shortening at its level
P = Ap * sigma_pm
sigma_cp = -P / Ac - P * zp * zp / Ic - M_unloaded * zp / Ic
tendon = dict(y=0.0, z=zp, Ap=Ap, eps_p0=keta.tendon_initial_strain(sigma_pm, sigma_cp, Ecm))

# by hand: transformed section (concrete + n Ap at the tendon) under M and the decompression force P0 at the tendon
n = Ep / Ecm
Ah = Ac + n * Ap
zh = n * Ap * zp / Ah
Ih = Ic + Ac * zh ** 2 + n * Ap * (zp - zh) ** 2
P0 = Ap * Ep * tendon["eps_p0"]

points = [(0.0, h / 2), (0.0, -h / 2), (0.0, zp)]     # top fibre, bottom fibre, tendon level
for label, M in (("unloaded", M_unloaded), ("in service", M_service)):
    state = keta.prestress_service_state(rings, [tendon], Ecm, sustained=(0.0, -M, 0.0), points=points)
    top, bottom, level = (s / 1e6 for s in state["point_stress"])
    print(f"{label}: top fibre {top:.2f} MPa, bottom fibre {bottom:.2f} MPa, concrete at tendon {level:.2f} MPa, "
          f"tendon {state['tendon_stress'][0] / 1e6:.1f} MPa")
    Mh = -M - P0 * zp + P0 * zh                       # moment about the transformed centroid
    by_hand = [-P0 / Ah + Mh * (z - zh) / Ih for _, z in points]
    print("  by hand:", ", ".join(f"{s / 1e6:.2f}" for s in by_hand), "MPa")
    assert all(abs(s - m) < 1e-9 * abs(m) for s, m in zip(state["point_stress"], by_hand))
    fmax, fmin = state["fibre_max"], state["fibre_min"]
    print(f"  extreme fibres: {fmax['stress'] / 1e6:.2f} MPa at z = {fmax['z']:.2f} m, "
          f"{fmin['stress'] / 1e6:.2f} MPa at z = {fmin['z']:.2f} m")
