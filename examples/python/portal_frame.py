"""A portal frame with load cases and combinations.

    cmake --build build --target _keta  # PYTHONPATH=build/python
    python examples/python/portal_frame.py
"""

import keta as k

SPAN, HEIGHT = 6.0, 3.5

m = k.Model()
a = m.add_node(0, 0, 0)
b = m.add_node(0, 0, HEIGHT)
c = m.add_node(SPAN, 0, HEIGHT)
d = m.add_node(SPAN, 0, 0)

steel = m.steel()
col_sec = m.add_section("HEB200")
beam_sec = m.add_section("IPE300")

col1 = m.add_beam(a, b, section=col_sec, material=steel)
beam = m.add_beam(b, c, section=beam_sec, material=steel)
col2 = m.add_beam(d, c, section=col_sec, material=steel)

# In 3D, a planar frame on true pins can still rotate rigidly about the
# axis through its two feet (a "trap-door" mechanism) -- keta rejects
# the singular system rather than inventing a result. Block rx at the feet:
# in-plane bending rotation (ry) stays free, so the frame behaves as pinned.
m.add_support(a, dx=True, dy=True, dz=True, rx=True)
m.add_support(d, dx=True, dy=True, dz=True, rx=True)

# Case 1: dead load (self-weight), case 2: live UDL, case 3: wind.
m.add_self_weight(case=1)
m.name_case(1, "Dead")
m.add_udl(beam, wz=-12_000, case=2)
m.name_case(2, "Live")
m.add_point_load(b, fx=8_000, case=3)
m.name_case(3, "Wind")

uls = m.add_combination([(1.35, 1), (1.5, 2)], name="ULS gravity")
uls_w = m.add_combination([(1.35, 1), (1.05, 2), (1.5, 3)], name="ULS wind lead")

r = m.analyze()

# Equilibrium by hand: the vertical reactions carry the factored self-weight
# (steel 7850 kg/m3, g = 9.80665 m/s2, areas of the catalogue) and live load,
# the horizontal ones the factored wind.
G, RHO, LIVE, WIND = 9.80665, 7850.0, 12_000.0 * SPAN, 8_000.0
weight = RHO * G * (2 * HEIGHT * k.section_props("HEB200")["A"] + SPAN * k.section_props("IPE300")["A"])
expected = {uls: (0.0, 1.35 * weight + 1.5 * LIVE), uls_w: (-1.5 * WIND, 1.35 * weight + 1.05 * LIVE)}
for combo_id, (fx, fz) in expected.items():
    sx, _, sz = r.combo(combo_id).sum_reactions()[:3]
    assert abs(sz - fz) < 1e-9 * fz, f"vertical equilibrium: {sz} vs {fz}"
    assert abs(sx - fx) < 1e-9 * fz, f"horizontal equilibrium: {sx} vs {fx}"

for combo_id, label in [(uls, "ULS gravity"), (uls_w, "ULS wind lead")]:
    res = r.combo(combo_id)
    print(f"\n=== {label} ===")
    print(f"max displacement: {res.max_displacement() * 1e3:.2f} mm")
    print(f"sum reactions   : {tuple(round(v / 1e3, 2) for v in res.sum_reactions())} kN")
    peaks = [max(abs(row[5]) for row in res.bar_forces(el)) for el in (col1, beam, col2)]
    for name, peak in zip(("column 1", "beam", "column 2"), peaks):
        print(f"  |My|max {name:9s}: {peak / 1e3:8.2f} kN.m")
