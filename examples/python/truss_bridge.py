"""A parametric Warren truss footbridge -- scripting as the second interface.

The whole geometry is generated from three parameters; change them and rerun.

    cmake --build build --target _keta  # PYTHONPATH=build/python
    python examples/python/truss_bridge.py
"""

import keta as k

N_BAYS = 8        # number of bottom-chord bays
BAY = 2.5         # m
RISE = 2.2        # truss height, m
LOAD = -15_000.0  # N per top node (deck + crowd)

m = k.Model()
steel = m.steel()
chord = m.add_section_tube(d=0.20, t=0.010)
diag = m.add_section_tube(d=0.14, t=0.008)

# Bottom and top chord nodes.
bottom = [m.add_node(i * BAY, 0, 0) for i in range(N_BAYS + 1)]
top = [m.add_node((i + 0.5) * BAY, 0, RISE) for i in range(N_BAYS)]

# A classic pin-jointed truss: every member axial-only. keta
# auto-eliminates the rotational DOFs of pure truss nodes (they carry no
# stiffness), so no fictitious rotational restraints are needed.
bars = []
for i in range(N_BAYS):
    bars.append(m.add_truss(bottom[i], bottom[i + 1], section=chord, material=steel))
for i in range(N_BAYS - 1):
    bars.append(m.add_truss(top[i], top[i + 1], section=chord, material=steel))
for i in range(N_BAYS):
    bars.append(m.add_truss(bottom[i], top[i], section=diag, material=steel))
    bars.append(m.add_truss(top[i], bottom[i + 1], section=diag, material=steel))

# Simply supported in-plane; the model is planar (X–Z), so block the
# out-of-plane translation dy at EVERY node -- this is a 3D solver.
m.pin(bottom[0])
m.add_support(bottom[-1], dy=True, dz=True)
for n in bottom[1:-1]:
    m.add_support(n, dy=True)
for n in top:
    m.add_support(n, dy=True)

m.add_self_weight(case=1)
m.name_case(1, "Dead")
for t in top:
    m.add_point_load(t, fz=LOAD, case=2)
m.name_case(2, "Crowd")
uls = m.add_combination([(1.35, 1), (1.5, 2)], name="ULS")

r = m.analyze()
res = r.combo(uls)

print(f"Warren truss: {N_BAYS} bays x {BAY} m = {N_BAYS * BAY} m span")
print(f"midspan deflection: {res.max_displacement() * 1e3:.1f} mm "
      f"(span/{N_BAYS * BAY / res.max_displacement():.0f})")

n_peak = max(
    (abs(row[1]), el) for el in bars for row in res.bar_forces(el)
)
print(f"|N|max = {n_peak[0] / 1e3:.1f} kN (bar {n_peak[1]})")
print(f"Sum Rz = {res.sum_reactions()[2] / 1e3:.1f} kN")

# The truss is statically determinate: under the crowd alone, the method of
# sections gives the middle bottom-chord bar. Cut bay MID and take moments
# about the top node above it: N RISE = R x - sum of the loads left of it.
crowd = r.case(2)
MID = N_BAYS // 2 - 1
pole_x = (MID + 0.5) * BAY
reaction = N_BAYS * abs(LOAD) / 2
moment = reaction * pole_x - sum(abs(LOAD) * (pole_x - (i + 0.5) * BAY) for i in range(MID))
n_chord = crowd.bar_forces(bars[MID])[0][1]
assert abs(n_chord - moment / RISE) < 1e-9 * moment / RISE, f"chord force {n_chord} vs {moment / RISE}"
assert abs(crowd.sum_reactions()[2] - N_BAYS * abs(LOAD)) < 1e-9 * N_BAYS * abs(LOAD)
