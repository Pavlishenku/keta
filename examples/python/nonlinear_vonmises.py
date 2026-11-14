"""keta example: snap-through of a shallow von Mises two-bar toggle (arc-length).

A shallow two-bar frame ("von Mises truss") pushed down at its apex is the
classic snap-through example: the load rises to a limit (fold) point, then the
apex snaps through to an inverted shape while the load *falls* — a branch a
load-controlled solver cannot follow (its tangent goes singular at the fold).
Cylindrical arc-length continuation traces the whole equilibrium path THROUGH
the fold, stamping the LDL^T inertia at every point and flagging the two
crossings where it flips (the fold, then the unloading-valley minimum).

The apex is restrained to the symmetric subspace (vertical only), so the discrete path lies on the hand-derived scalar
curve P(w) and its fold is the closed-form limit load. The reference load is
set to that closed-form fold, so the trace's load peak sits at lambda ~ 1.

Run after building the extension:
    cmake --build build --target _keta   (then PYTHONPATH=build/python)
    python examples/python/nonlinear_vonmises.py
"""

import math

import keta as k

# Frame geometry (SI: m, N) — a half-span C, rise H shallow arch.
C = 4.0            # half-span (m)
H = 1.0            # apex rise (m)
E = 210e9          # steel Young's modulus (Pa)
SX = 1.0e-2        # section area A (m^2)
IY = 1.0e-5        # second moment of area (m^4)
IX = 2.0e-5        # torsion constant (m^4)


def apex_load(w):
    """Hand-derived apex load P(w) for downward apex descent w (the closed form
    the symmetric discrete path lies on)."""
    l0 = math.hypot(C, H)
    hn = H - w
    ln = math.hypot(C, hn)
    delta = math.atan(H / C) - math.atan2(hn, C)
    return 2.0 * (
        (E * SX / l0) * (ln - l0) * (-hn / ln)
        + (12.0 * E * IY / l0) * math.sin(delta) * math.cos(delta) * C / (ln * ln)
    )


def fold_load():
    """The limit (fold) load: maximise P(w) on [0, H] by ternary search."""
    lo, hi = 0.0, H
    for _ in range(300):
        m1 = lo + (hi - lo) / 3.0
        m2 = hi - (hi - lo) / 3.0
        if apex_load(m1) < apex_load(m2):
            lo = m1
        else:
            hi = m2
    w = 0.5 * (lo + hi)
    return w, apex_load(w)


w_lim, p_lim = fold_load()

# Build the frame: two beams from the feet up to the apex.
m = k.Model()
steel = m.steel()                                   # E = 210 GPa
sec = m.add_section_props(sx=SX, iy=IY, iz=IY, ix=IX)
n1 = m.add_node(-C, 0.0, 0.0)
n2 = m.add_node(C, 0.0, 0.0)
apex = m.add_node(0.0, 0.0, H)
m.add_beam(n1, apex, section=sec, material=steel)
m.add_beam(n2, apex, section=sec, material=steel)

# Clamp both feet; restrain the apex to the symmetric subspace (vertical only).
m.fix(n1)
m.fix(n2)
m.add_support(apex, dx=True, dy=True, rx=True, ry=True, rz=True)

# Reference apex load = the closed-form fold, so the load peak is at lambda ~ 1.
m.add_point_load(apex, fz=-p_lim, case=1)

problems = m.validate()
assert not problems, problems

# Trace the equilibrium path through the fold. arc_length is |delta w| on the
# 1-DOF symmetric subspace; stop on the far rising branch past the snap-through.
trace = k.Model.arc_length(
    m, case=1, arc_length=0.01, max_steps=600, stop_lambda=1.3
)
print(trace)

lambdas = trace.lambdas()
crossings = trace.crossings()

# The fold is a LOCAL load peak on the first branch: past the snap-through the
# far rising branch climbs higher (that is where stop_lambda halts the trace),
# so the fold is the maximum lambda BEFORE the unloading-valley crossing.
valley = crossings[1] if len(crossings) >= 2 else len(lambdas)
peak_i = max(range(valley), key=lambda i: lambdas[i])

# lambda vs apex-descent w = -uz(apex), sampled through the fold (plus the exact
# fold and crossing rows).
print("\n   pt      w (m)     lambda    n_neg   note")
print("  ---------------------------------------------")
marks = set(range(0, trace.n_points, 12)) | set(crossings) | {peak_i}
for i in sorted(marks):
    w = -trace.node_disp(i, apex)[2]
    note = ""
    if i == peak_i:
        note = "<- fold (limit load)"
    elif i in crossings:
        note = "<- inertia crossing"
    print(f"  {i:4d}   {w:8.4f}   {lambdas[i]:8.4f}   {trace.n_neg(i):3d}    {note}")

# The fold: the load peak on the first branch (a local maximum in lambda).
w_peak = -trace.node_disp(peak_i, apex)[2]
print(f"\nload peak (fold)    : lambda = {lambdas[peak_i]:.4f} at w = {w_peak:.4f} m")
print(f"closed-form fold    : lambda = 1.0000 at w = {w_lim:.4f} m  (P_lim = {p_lim/1e3:.1f} kN)")
print(f"minimum lambda      : {min(lambdas):.4f}   (unloading valley, load reversed)")
print(f"inertia crossings   : {crossings}")
for c in crossings:
    wc = -trace.node_disp(c, apex)[2]
    print(f"   at point {c:3d}: w = {wc:.4f} m, lambda = {lambdas[c]:.4f}, "
          f"n_neg {trace.n_neg(c - 1)} -> {trace.n_neg(c)}")

# A genuine snap-through fold: load rises, peaks (interior), falls below zero,
# recovers on the inverted branch — with exactly two inertia crossings.
assert 0 < peak_i < valley - 1, "fold peak must be interior on the first branch"
assert abs(lambdas[peak_i] - 1.0) < 0.05, "peak lambda ~ 1 at the closed-form fold"
assert min(lambdas) < 0.0, "the unloading valley must reverse the load"
assert len(crossings) == 2, f"a snap-through fold has two crossings, got {crossings}"
print("\nOK: the arc-length trace folds through the limit point (rise, peak, "
      "fall, recover) with two inertia crossings.")
