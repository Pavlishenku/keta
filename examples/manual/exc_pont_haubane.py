# COMPLETE EXAMPLE 2 - Cable-stayed bridge with a central pylon (spatial structure):
#   box girder deck of 2 x 60 m supported on the abutments and on the pylon; pylon of 50 m
#   (20 m below the deck); 2 x 4 stays in a fan arrangement in the axial plane;
#   tuning of the stays under permanent loads, mode shapes and live load.
import keta

m = keta.Model()
concrete, steel = m.add_material(e=36e9, nu=0.2, rho=2500.0), m.steel()
box = m.add_section_props(sx=6.0, iy=4.0, iz=60.0, ix=8.0, name="box")
shaft = m.add_section_props(sx=8.0, iy=10.0, iz=10.0, ix=15.0, name="pylon")
stay_sec = m.add_section_props(sx=6.0e-3, iy=1e-9, iz=1e-9, ix=1e-9, name="stay")

xs = [-60.0 + 6.0 * i for i in range(21)]
deck = {x: m.add_node(x, 0.0, 0.0) for x in xs}
deck_bars = [m.add_beam(deck[xs[i]], deck[xs[i + 1]], box, concrete) for i in range(20)]
pylon = [m.add_node(0.0, 0.0, z) for z in (-20.0, -10.0, 0.0, 10.0, 20.0, 30.0)]
pylon_bars = [m.add_beam(pylon[i], pylon[i + 1], shaft, concrete) for i in range(5)]
head = pylon[-1]
anchors = [-48.0, -36.0, -24.0, -12.0, 12.0, 24.0, 36.0, 48.0]
stays = [m.add_truss(head, deck[x], stay_sec, steel) for x in anchors]

m.fix(pylon[0])                                                     # clamped at the pylon base
m.add_support(deck[-60.0], dy=True, dz=True, rx=True)               # abutments: sliding supports
m.add_support(deck[60.0], dy=True, dz=True, rx=True)
m.add_rigid_link(pylon[2], [deck[0.0]], dx=True, dy=True, dz=True, rx=True, ry=False, rz=False)  # support on the pylon

m.add_self_weight(case=1)
for b in deck_bars:
    m.add_udl(b, wz=-60e3, case=1)                                  # superimposed dead loads
unit_cases = [m.add_unit_prestress_case(h, P=1.0) for h in stays]   # cases 2 to 9
m.name_case(1, "permanent")

# Tuning: zero deflection under permanent loads at the eight anchorage points
targets = [dict(kind="displacement", id=deck[x], comp=2, value=0.0) for x in anchors]
t = m.tune_linear(base=[1], unknowns=unit_cases, targets=targets)
print("tuning tensions (unit cases) [MN]:", [round(v / 1e6, 4) for v in t["x"]])
tuned = m.add_combination([(1.0, 1)] + [(v, c) for v, c in zip(t["x"], unit_cases)], name="tuned permanent")

# Live load: 3 lanes (LM1 distributed 9 + 2 x 2.5 kN/m2 x 3 m) on the right span
m.add_point_load(deck[0.0], fz=0.0, case=20)
for b in deck_bars[10:]:
    m.add_udl(b, wz=-(9e3 + 2 * 2.5e3) * 3.0, case=20)
m.name_case(20, "live load right span")
r = m.analyze(stations=5)
perm = r.combo(tuned)
fmax = max(abs(perm.displacement(deck[x])[2]) for x in anchors)
assert fmax < 1e-9
print("tuned permanent: zero deflections at the anchorages (|uz| < 1e-9 m)")
print("forces in the stays [kN]:", [round(perm.bar_forces(h)[0][1] / 1e3, 1) for h in stays])
print(f"pylon: N at the base {perm.bar_forces(pylon_bars[0])[0][1] / 1e3:.1f} kN, "
      f"head displacement dx = {round(perm.displacement(head)[0] * 1e3, 3) + 0.0:.3f} mm")
q = r.case(20)
print(f"live load: maximum deflection {min(q.displacement(n)[2] for n in deck.values()) * 1e3:.2f} mm, "
      f"dx of the pylon head {q.displacement(head)[0] * 1e3:.2f} mm, "
      f"max tension variation {max(abs(q.bar_forces(h)[0][1]) for h in stays) / 1e3:.1f} kN")

mo = m.modal(n_modes=5)
for k in range(5):
    f = mo.frequencies[k]
    fr = mo.mass_ratio()[k]
    print(f"mode {k + 1}: {f:6.3f} Hz, effective masses X {fr[0] * 100:5.1f} %  Y {fr[1] * 100:5.1f} %  "
          f"Z {fr[2] * 100:5.1f} %")
