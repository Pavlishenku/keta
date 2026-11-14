"""keta example: modal (free-vibration) analysis of a simply-supported steel beam.

A slender 6 m steel beam, pinned at one end and on a roller at the other, is a
textbook simply-supported bending problem. Its natural frequencies are

    f_n = (n*pi)^2 / (2*pi*L^2) * sqrt(E*I / (rho*A)),   n = 1, 2, 3, ...

The fundamental bends about the weak principal axis (the smaller second moment
of area), so f_1 uses I_min. This script runs the modal solver, prints the first
three frequencies and checks f_1 against the closed form (< 1% relative error).

Modal analysis needs mass: a material must carry a density. The steel() preset
already sets rho = 7850 kg/m^3, so no extra setup is required.

Run after building the extension:
    cmake --build build --target _keta   (then PYTHONPATH=build/python)
    python examples/python/modal_beam.py
"""

import math

import keta as k

L = 6.0            # beam length (m)
N_ELEMS = 12       # mesh density (consistent mass converges quickly; even count
                   # so a node lands exactly at midspan)
B, H = 0.05, 0.15  # rectangular section: b along local y, h along local z (m)

m = k.Model()
steel = m.steel()               # E = 210 GPa, nu = 0.3, rho = 7850 kg/m^3
rect = m.add_section_rect(B, H)

# The beam runs along +X, meshed into N_ELEMS beam elements.
nodes = [m.add_node(L * i / N_ELEMS, 0.0, 0.0) for i in range(N_ELEMS + 1)]
for a, b in zip(nodes[:-1], nodes[1:]):
    m.add_beam(a, b, section=rect, material=steel)

# Simply-supported end conditions. The base is a pin (all translations held) and
# also holds twist about the axis (rx) so the beam has no rigid-body torsional
# mode. The far end is a roller: transverse translations held, but it slides
# axially. Bending rotations (ry, rz) stay free at both ends -> pinned-pinned
# flexure, a clean half sine.
m.add_support(nodes[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(nodes[-1], dy=True, dz=True)

problems = m.validate()
assert not problems, problems

# Modal analysis: the six lowest natural frequencies and their mode shapes.
modal = m.modal(n_modes=6)
print(modal)

print("\nfirst three natural frequencies:")
for i, f in enumerate(modal.frequencies[:3], start=1):
    print(f"  mode {i}: f = {f:8.3f} Hz   (T = {1.0 / f:6.4f} s)")

# Closed form for the fundamental. The rectangle's second moments of area are
# I_y = b*h^3/12 (bending in x-z) and I_z = h*b^3/12 (bending in x-y); the beam
# vibrates first about the weak axis, so f_1 uses the smaller one.
E, RHO = 210e9, 7850.0
area = B * H
i_y = B * H**3 / 12.0
i_z = H * B**3 / 12.0
i_min = min(i_y, i_z)
f1_closed = (math.pi**2) / (2.0 * math.pi * L**2) * math.sqrt(E * i_min / (RHO * area))

f1_fem = modal.frequencies[0]
rel_err = abs(f1_fem - f1_closed) / f1_closed
print(f"\nf_1 (FEM)        : {f1_fem:9.4f} Hz")
print(f"f_1 (closed form): {f1_closed:9.4f} Hz")
print(f"relative error   : {rel_err * 100:9.4f} %")
assert rel_err < 0.01, f"fundamental frequency off by {rel_err * 100:.3f}%"

# Fundamental mode shape: the unit-crest half sine peaks at midspan, deflecting
# in the weak-axis direction (global Y here). mode_displacement returns the full
# (dx, dy, dz, rx, ry, rz) tuple for a node in a given (0-based) mode.
mid = nodes[N_ELEMS // 2]
dx, dy, dz, rx, ry, rz = modal.mode_displacement(0, mid)
print("\nmode 1 midspan displacement (unit-crest):")
print(f"  dy = {dy:+.4f}   dz = {dz:+.4f}   (dy is the +1 crest, weak-axis bending)")
assert abs(dy) > 0.99, "the fundamental crest should sit at midspan in global Y"
