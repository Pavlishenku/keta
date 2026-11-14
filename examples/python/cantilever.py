"""keta quickstart: a steel cantilever with a tip load.

Run after building the extension:
    cmake --build build --target _keta  # PYTHONPATH=build/python
    python examples/python/cantilever.py
"""

import keta as k

L = 4.0  # m
P = -10_000.0  # N (downward)

m = k.Model()
a = m.add_node(0, 0, 0)
b = m.add_node(L, 0, 0)

steel = m.steel()
ipe = m.add_section("IPE300")
beam = m.add_beam(a, b, section=ipe, material=steel)

m.fix(a)
m.add_point_load(b, fz=P, case=1)

problems = m.validate()
assert not problems, problems

r = m.analyze()
case = r.case(1)

ux, uy, uz, rx, ry, rz = case.displacement(b)
print(f"tip deflection : {uz * 1e3:8.3f} mm")
print(f"tip rotation   : {ry:8.6f} rad")

# Compare with the Euler-Bernoulli closed form (shear makes it slightly softer).
props = k.section_props("IPE300")
E, Iy = 210e9, props["Iy"]
print(f"EB closed form : {P * L**3 / (3 * E * Iy) * 1e3:8.3f} mm")

# The tip rotation of a Timoshenko cantilever under a tip load has no shear
# part: |theta| = |P| L^2 / (2 E I) exactly. The deflection adds the shear
# term to the Euler-Bernoulli value, so it is larger, by a few percent here.
theta = abs(P) * L**2 / (2 * E * Iy)
assert abs(abs(ry) - theta) < 1e-9 * theta, f"tip rotation {ry} vs {theta}"
eb = abs(P) * L**3 / (3 * E * Iy)
assert eb < abs(uz) < 1.05 * eb, f"tip deflection {uz} vs Euler-Bernoulli {eb}"

print("\nsupport reactions:")
for node, (fx, fy, fz, mx, my, mz) in case.reactions():
    print(f"  node {node}: Fz = {fz / 1e3:.2f} kN, My = {my / 1e3:.2f} kN.m")
    # statics of the cantilever: the clamp carries |P| and |P| L
    assert abs(fz + P) < 1e-9 * abs(P) and abs(abs(my) - abs(P) * L) < 1e-9 * abs(P) * L

print("\nbending moment along the beam:")
for xi, n, vy, vz, mt, my, mz in case.bar_forces(beam):
    bar_txt = "#" * int(abs(my) / abs(P * L) * 40)
    print(f"  x/L = {xi:4.2f}  My = {my / 1e3:9.2f} kN.m  {bar_txt}")
