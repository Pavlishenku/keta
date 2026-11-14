"""keta example: linear (eigenvalue) buckling of a pinned-pinned steel column.

A slender 4 m column under a unit axial compression. Linear buckling returns
the critical load factors λ that scale the reference (unit) load case to the
buckling loads, so here λ₁ is the Euler critical force itself and must match
the closed form π²·E·I_min / L² (buckling about the weak principal axis).

Run after building the extension:
    cmake --build build --target _keta   (then PYTHONPATH=build/python)
    python examples/python/buckling_column.py
"""

import math

import keta as k

L = 4.0            # column length (m)
N_ELEMS = 10       # mesh density (buckling needs a few elements to converge)
B, H = 0.06, 0.12  # rectangular section: b along local y, h along local z (m)
P_REF = 1.0        # unit reference compression (N) → λ₁ is P_cr in newtons

m = k.Model()
steel = m.steel()
rect = m.add_section_rect(B, H)

# The column runs along +X, meshed into N_ELEMS beam elements.
nodes = [m.add_node(L * i / N_ELEMS, 0.0, 0.0) for i in range(N_ELEMS + 1)]
for a, b in zip(nodes[:-1], nodes[1:]):
    m.add_beam(a, b, section=rect, material=steel)

# Pinned-pinned end conditions: the base is a full pin (all translations held)
# and also holds twist about the axis (rx) — otherwise the column, free to
# rotate about its own axis at every node, has a rigid-body torsional mode. The
# top is a lateral pin that still slides axially so the reference load actually
# compresses the member. Bending rotations (ry, rz) stay free at both ends, so
# the flexural mode is a clean half sine.
m.add_support(nodes[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(nodes[-1], dy=True, dz=True)

# Unit axial compression at the top, directed toward the base.
m.add_point_load(nodes[-1], fx=-P_REF, case=1)

problems = m.validate()
assert not problems, problems

# Linear buckling on load case 1 (the reference compression).
buck = m.buckling(case=1, n_modes=4)
print(buck)

lam1 = buck.factors[0]

# Euler closed form: a pinned-pinned column buckles about the weak (smaller-I)
# principal axis. For the rectangle, I about local y is B·H³/12 and about
# local z is H·B³/12; the smaller governs.
E = 210e9
i_min = min(B * H**3 / 12.0, H * B**3 / 12.0)
euler = math.pi**2 * E * i_min / L**2

rel_err = abs(lam1 - euler) / euler
print(f"\nlambda_1 (FEM)   : {lam1:12.1f} N")
print(f"Euler pi^2EI/L^2 : {euler:12.1f} N")
print(f"relative error   : {rel_err * 100:9.4f} %")
assert rel_err < 0.01, f"buckling factor off by {rel_err * 100:.3f}%"

# First mode's midspan crest — the unit-crest half sine peaks at the middle
# node, deflecting in the weak-axis (global Y) direction.
mid = nodes[N_ELEMS // 2]
dx, dy, dz, rx, ry, rz = buck.mode_displacement(0, mid)
print("\nmode 1 midspan displacement (unit-crest):")
print(f"  dy = {dy:+.4f}   dz = {dz:+.4f}")

print("\nall critical load factors:")
for i, lam in enumerate(buck.factors, start=1):
    print(f"  mode {i}: lambda = {lam:12.1f} N")
