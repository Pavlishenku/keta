"""keta example: response-spectrum (spectral) seismic analysis of a shear building.

A two-storey shear building is the textbook two-DOF earthquake problem. Each
storey is idealized as a spring of stiffness k = E*A/L acting along the ground-
motion direction (+X here), with the floor masses lumped at the two levels. In
keta this is a clean axial chain: a massless bar per storey (the storey
spring) with a lumped point mass at each floor, every DOF but the axial one held.
The stiffness / mass matrices are then

    K = [[k1 + k2, -k2], [-k2, k2]]        M = diag(m1, m2)

For a unit ground acceleration in X, mode i (M-normalized, phi_i^T M phi_i = 1)
has participation Gamma_i = phi_i^T M r (r = [1, 1]), spectral displacement
Sd_i = Sa(T_i) / w_i^2, and a signed peak floor field u_i = Gamma_i * Sd_i * phi_i.
The floors do not peak together, so the modal peaks are combined by SRSS,
u = sqrt(sum_i u_i^2) per floor -- never added.

The script runs the response-spectrum solver against a flat design spectrum and
checks:

  * the SRSS roof displacement matches a fully independent 2x2 hand solution;
  * the effective modal masses Gamma_i^2 sum to the total directional mass
    m1 + m2 (the completeness relation, with every mode retained).

Run after building the extension:
    cmake --build build --target _keta   (then PYTHONPATH=build/python)
    python examples/python/response_spectrum_building.py
"""

import math

import keta as k

E = 210e9          # Young's modulus (Pa)
L = 3.0            # storey height (m)
K_STOREY = 20e6    # target storey stiffness k = E*A/L (N/m)
MASS = 10_000.0    # floor mass, both floors equal (kg)
SA = 2.0           # flat design pseudo-acceleration (m/s^2)
ZETA = 0.05        # modal damping (5 %, the usual seismic value)

area = K_STOREY * L / E     # pick A so E*A/L is exactly K_STOREY

m = k.Model()
mat = m.add_material(E, nu=0.3)                       # massless storey springs
sec = m.add_section_props(sx=area, iy=1e-5, iz=1e-5, ix=1e-5)

# Base (fixed) and two floors, stacked along +X (the ground-motion direction).
base = m.add_node(0.0, 0.0, 0.0)
floor1 = m.add_node(L, 0.0, 0.0)
floor2 = m.add_node(2.0 * L, 0.0, 0.0)
m.add_beam(base, floor1, section=sec, material=mat)   # storey 1 spring
m.add_beam(floor1, floor2, section=sec, material=mat)  # storey 2 spring

# Shear-building kinematics: floors translate only in X; the base is clamped.
m.fix(base)
for n in (floor1, floor2):
    m.add_support(n, dy=True, dz=True, rx=True, ry=True, rz=True)
m.add_point_mass(floor1, MASS)
m.add_point_mass(floor2, MASS)

problems = m.validate()
assert not problems, problems

# A flat spectrum: two equal ordinates make Sa(T) = SA at every period.
spec = m.response_spectrum(
    spectrum=[(0.01, SA), (10.0, SA)],
    dir=0,           # ground motion along global X
    zeta=ZETA,
    rule="srss",
)
print(spec)

# --- Modal / effective-mass table ------------------------------------------
freqs = spec.frequencies()
part = spec.participation()
eff = spec.effective_mass()
total = spec.total_mass()
print("\nmode   f (Hz)   participation   eff. mass (kg)   ratio")
for i, (f, g, e) in enumerate(zip(freqs, part, eff), start=1):
    print(f"  {i}   {f:7.3f}   {g:12.4f}   {e:12.1f}   {e / total:6.3f}")
print(f"                                sum = {sum(eff):9.1f}   (total = {total:.1f})")
assert abs(sum(eff) - total) / total < 1e-9, "effective masses must sum to the total"
assert abs(total - 2.0 * MASS) / (2.0 * MASS) < 1e-12, "total directional mass = m1 + m2"

# --- SRSS roof displacement vs an independent 2x2 hand solution -------------
roof_fem = spec.combined_disp(floor2)[0]     # global-X roof displacement
floor1_fem = spec.combined_disp(floor1)[0]


def hand_solution():
    """Independent SRSS roof/floor displacement of the 2-DOF chain."""
    k1 = k2 = K_STOREY
    # Eigenvalues of M^-1 K (M = MASS * I) via the 2x2 closed form.
    a, b, d = (k1 + k2) / MASS, -k2 / MASS, k2 / MASS
    tr, det = a + d, a * d - b * b
    root = math.sqrt(tr * tr - 4.0 * det)
    lams = [(tr - root) / 2.0, (tr + root) / 2.0]
    u = [0.0, 0.0]
    for lam in lams:
        w = math.sqrt(lam)
        x, y = 1.0, (lam - a) / b            # mode shape (a - lam) x + b y = 0
        norm = math.sqrt(MASS * (x * x + y * y))  # M-normalize
        x, y = x / norm, y / norm
        gamma = MASS * (x + y)               # phi^T M r, r = [1, 1]
        sd = SA / (w * w)
        u[0] += (gamma * sd * x) ** 2        # SRSS accumulation
        u[1] += (gamma * sd * y) ** 2
    return math.sqrt(u[0]), math.sqrt(u[1])


floor1_hand, roof_hand = hand_solution()
print(f"\nSRSS floor-1 drift : {floor1_fem * 1e3:8.4f} mm   (hand {floor1_hand * 1e3:.4f} mm)")
print(f"SRSS roof drift    : {roof_fem * 1e3:8.4f} mm   (hand {roof_hand * 1e3:.4f} mm)")
print(f"base shear (SRSS)  : {spec.base_shear() / 1e3:8.4f} kN")
rel = abs(roof_fem - roof_hand) / roof_hand
print(f"roof relative error: {rel * 100:.3e} %")
assert rel < 1e-9, f"SRSS roof displacement off by {rel:.2e}"
assert abs(floor1_fem - floor1_hand) / floor1_hand < 1e-9

print("\nOK: SRSS roof displacement matches the hand-derived 2-DOF spectrum, and")
print("    the effective modal masses sum to the total storey mass.")
