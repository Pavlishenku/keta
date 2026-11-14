"""keta example: a cable-stayed cantilever with an exact elastic-catenary stay.

A vertical steel mast and a horizontal deck cantilever share a clamped base. A
prestressed stay cable is slung from the mast top down to the deck tip, and a
downward point load hangs at the tip. The stay is an exact elastic-catenary
element: it sags under its own weight and stiffens as the tip pulls on it, so the
structure is geometrically nonlinear and is solved with the dedicated cable
driver. Its unstrained length is derived from a prescribed horizontal tension H
(the prestress / form-finding route), resolved once at the undeformed geometry.

The script compares the deck-tip deflection WITH the stay against the cable-less
linear cantilever, and prints the stay's horizontal tension H, its maximum
tension T_max and its mid-span sag.

Run after building the extension:
    cmake --build build --target _keta   (then PYTHONPATH=build/python)
    python examples/python/cable_stayed_mast.py
"""

import math

import keta as k

MAST_H = 6.0        # mast height (m)
DECK_L = 8.0        # deck cantilever length (m)
B = 0.15            # solid square section b x b (m) for the mast + deck
P_TIP = -10_000.0   # downward point load at the deck tip (N)
EA = 1.5e8          # stay axial rigidity EA (N)
W_CABLE = 800.0     # stay self-weight per unstrained length (N/m)
H_STAY = 6_000.0    # prescribed horizontal tension of the stay (N) -> the prestress

m = k.Model()
steel = m.steel()
sq = m.add_section_rect(B, B)

# base -> clamped, mast top, deck tip.
base = m.add_node(0.0, 0.0, 0.0)
top = m.add_node(0.0, 0.0, MAST_H)
tip = m.add_node(DECK_L, 0.0, 0.0)
m.add_beam(base, top, section=sq, material=steel)   # the mast
m.add_beam(base, tip, section=sq, material=steel)   # the deck cantilever
m.fix(base)
m.add_point_load(tip, fz=P_TIP, case=1)

problems = m.validate()
assert not problems, problems

# Cable-less reference: the plain linear cantilever tip deflection (the mast is
# unloaded and only shares the clamped base, so it contributes nothing here).
# `analyze` is run BEFORE the cable is added (the linear driver rejects a cabled
# model).
lin = m.analyze().case(1)
tip_lin = lin.displacement(tip)[2]   # global-Z deflection (m), negative = down

# Add the prestressed elastic-catenary stay from the mast top to the deck tip.
# The horizontal-tension route derives the unstrained length so H = H_STAY at the
# undeformed chord; self-weight w then sags it.
stay = m.add_cable(top, tip, ea=EA, w=W_CABLE, h=H_STAY)

res = m.cable_static(case=1)
print(res)
tip_cab = res.displacement(tip)[2]   # global-Z deflection WITH the stay (m)

# Stay end forces and the sampled catenary polyline.
h, v1, v2, t_max, l0 = res.cable(stay)
shape = res.cable_shape(stay)

# Mid-span sag: the largest drop of the sampled centre-line below the straight
# chord between the two (deformed) anchor points.
a, z_a = shape[0], shape[0][2]
b, z_b = shape[-1], shape[-1][2]
n = len(shape) - 1
sag = max(
    (z_a + (z_b - z_a) * (i / n)) - p[2]   # chord_z - point_z (down = positive)
    for i, p in enumerate(shape)
)

chord = math.hypot(DECK_L, MAST_H)
print(f"\nstay unstrained length L0 : {l0:8.3f} m   (undeformed chord approx {chord:.3f} m)")
print(f"prescribed H (undeformed) : {H_STAY / 1e3:8.2f} kN   (sets L0 at the reference geometry)")
print(f"converged H (deformed)    : {h / 1e3:8.2f} kN   (higher: the sagging tip re-tensions the stay)")
print(f"end tensions V1 / V2      : {v1 / 1e3:8.2f} / {v2 / 1e3:.2f} kN")
print(f"maximum tension T_max     : {t_max / 1e3:8.2f} kN")
print(f"mid-span sag              : {sag * 1e3:8.1f} mm")
print(f"\ndeck-tip deflection (no stay) : {tip_lin * 1e3:8.2f} mm")
print(f"deck-tip deflection (stayed)  : {tip_cab * 1e3:8.2f} mm   ({res.iterations} iterations)")
reduction = (abs(tip_lin) - abs(tip_cab)) / abs(tip_lin)
print(f"stay reduces the tip sag by   : {reduction * 100:8.1f} %")

# The stay carries a real tension state (T_max >= H > 0). The prescribed H fixes
# the unstrained length at the undeformed chord; the converged H differs once the
# deck deflects and re-tensions the stay -- so L0 stays near the undeformed chord.
assert t_max >= h > 0.0, f"tension state invalid: H={h}, T_max={t_max}"
assert abs(l0 - chord) < 0.1 * chord, f"L0 should track the undeformed chord: {l0} vs {chord}"
# The stay is a working stay: it lifts the tip, so the deck sags less than the
# bare cantilever.
assert abs(tip_cab) < abs(tip_lin), "the stay must reduce the tip deflection"
# The self-weight cable actually sags (a taut string would not).
assert sag > 1e-3, f"the stay should sag under self-weight: {sag} m"
print("\nOK: the prestressed elastic-catenary stay sags under self-weight and lifts the deck tip.")
