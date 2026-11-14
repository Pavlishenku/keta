"""keta example: second-order (P-Delta) sway amplification of a cantilever column.

A slender steel column, fixed at the base, carries a heavy axial compression at
its free top together with a small lateral push. The axial force adds a negative
geometric stiffness that softens the column, amplifying the first-order sway by
the classic beam-column factor 1/(1 - P/P_cr). This script compares the linear
and the second-order tip drift and checks the amplification against that closed
form (and against the engine's own linear buckling load).

Run after building the extension:
    cmake --build build --target _keta   (then PYTHONPATH=build/python)
    python examples/python/second_order_column.py
"""

import math

import keta as k

H = 5.0            # column height (m)
N_ELEMS = 10       # mesh density (a few elements converge the P-Delta curve)
B = 0.1            # solid square section b x b (m): Iy = Iz, sway is isotropic
P = 100_000.0      # axial compression at the top (N) -> ~0.58 * P_cr
Q = 3_000.0        # small lateral push at the top (N)

m = k.Model()
steel = m.steel()
sq = m.add_section_rect(B, B)

# The column runs up +Z, meshed into N_ELEMS beam elements, fixed at the base.
nodes = [m.add_node(0.0, 0.0, H * i / N_ELEMS) for i in range(N_ELEMS + 1)]
for a, b in zip(nodes[:-1], nodes[1:]):
    m.add_beam(a, b, section=sq, material=steel)
m.fix(nodes[0])
top = nodes[-1]

# One load case: axial compression (down, -Z) plus a small lateral push (+X) at
# the free top. The axial load alone causes no drift; the lateral load sets the
# first-order sway that the axial force then amplifies.
m.add_point_load(top, fz=-P, fx=Q, case=1)

problems = m.validate()
assert not problems, problems

# First-order (linear) tip drift — the global-X translation of the free top.
lin = m.analyze().case(1)
drift_lin = abs(lin.displacement(top)[0])

# Second-order (P-Delta) tip drift, same load case.
so = m.second_order(case=1)
drift_so = abs(so.displacement(top)[0])
print(so)

amp_fem = drift_so / drift_lin

# Closed-form Euler critical load of a cantilever (effective length 2H) and the
# beam-column amplification factor.
E = 210e9
I = B**4 / 12.0
Pcr = math.pi**2 * E * I / (4.0 * H**2)
ratio = P / Pcr
amp_theory = 1.0 / (1.0 - ratio)

# Cross-check P_cr against the engine's own linear buckling factor: lambda_1
# scales the reference axial load to the critical load.
lam1 = m.buckling(case=1, n_modes=1).factors[0]
Pcr_fem = lam1 * P

print(f"\naxial P             : {P / 1e3:8.1f} kN")
print(f"Euler P_cr (form)   : {Pcr / 1e3:8.1f} kN   (P/P_cr = {ratio:.3f})")
print(f"P_cr (FEM buckling) : {Pcr_fem / 1e3:8.1f} kN   (lambda_1 = {lam1:.4f})")
print(f"\ntip drift linear    : {drift_lin * 1e3:8.2f} mm")
print(f"tip drift 2nd-order : {drift_so * 1e3:8.2f} mm   ({so.iterations} iterations)")
print(f"\namplification (FEM) : {amp_fem:.4f}")
print(f"1/(1 - P/P_cr)      : {amp_theory:.4f}")
rel = abs(amp_fem - amp_theory) / amp_theory
print(f"relative difference : {rel * 100:.2f} %")

assert rel < 0.05, f"amplification off by {rel * 100:.2f}%"
print("\nOK: second-order sway matches 1/(1 - P/P_cr) within a few percent.")
