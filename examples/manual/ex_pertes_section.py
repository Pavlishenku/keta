# Rectangular beam 0.40 x 0.80 m in C40/50 (Ecm = 35 GPa), two layers of strands pretensioned to
# 1400 MPa on the bed: state at release (elastic shortening), time-dependent losses at
# 90, 365 and 20,000 days, then resisting moment before and after losses.
import math
import keta

b, h, fck, Ecm, Ep = 0.40, 0.80, 40e6, 35e9, 195e9
fcd = 0.85 * fck / 1.5
rings = [[(-b / 2, -h / 2), (b / 2, -h / 2), (b / 2, h / 2), (-b / 2, h / 2)]]
bed = 1400e6                                              # stress on the bed: eps_p0 = bed / Ep
layers = [dict(y=0.0, z=-0.32, Ap=900e-6, eps_p0=bed / Ep), dict(y=0.0, z=-0.25, Ap=600e-6, eps_p0=bed / Ep)]

state = keta.prestress_service_state(rings, layers, Ecm)
for layer, sp, sc in zip(layers, state["tendon_stress"], state["concrete_stress"]):
    print(f"z = {layer['z']:.2f} m: strands {sp / 1e6:.1f} MPa (elastic loss {(bed - sp) / 1e6:.1f}), "
          f"concrete {sc / 1e6:.2f} MPa")
# by hand for a single layer of the total area at the centroid of the strands: sigma = bed / (1 + n Ap (1/Ac + e²/Ic))
Ap, e = 1500e-6, (900 * 0.32 + 600 * 0.25) / 1500
print(f"  equivalent single layer: {bed / (1 + Ep / Ecm * Ap * (1 / (b * h) + e * e / (b * h ** 3 / 12))) / 1e6:.1f} MPa")

history = keta.section_loss_history(rings, layers, [90.0, 365.0, 20000.0], fck, rh=70.0, h0=0.267, t0=28.0, Ecm=Ecm)
for age, losses in history:
    print(f"{age:7.0f} days:", ", ".join(
        f"{p['total'] / 1e6:6.1f} MPa (creep {p['creep'] / 1e6:.1f}, shrinkage {p['shrinkage'] / 1e6:.1f}, "
        f"relaxation {p['relaxation'] / 1e6:.1f})" for p in losses))

final = history[-1][1]
after = [dict(layer, eps_p0=p["eps_p0_after"]) for layer, p in zip(layers, final)]
for label, tendons in (("before losses", layers), ("after losses", after)):
    r = keta.prestressed_section_resistance(rings, fck, fcd, tendons, theta=math.pi)
    print(f"MRd {label}: {r['M_R'] / 1e3:.1f} kN.m")
print(f"post-tensioning, 4 tendons stressed one after the other: mean share {keta.sequential_tensioning_factor(4):.3f}")
