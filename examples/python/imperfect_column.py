"""keta example: mode-based imperfections and nodal point masses.

Two just-added engine features, exercised end to end against closed forms.

1. Imperfect column (EN 1993-1-1 section 5.3.2). A pinned-pinned steel column is
   buckled, and its lowest buckling mode is re-injected as an initial geometric
   imperfection scaled to an erection-tolerance bow delta0 = L/1000. A
   second-order (P-Delta) solve at a working compression P = 0.8 * P_cr then
   amplifies that bow. The load-induced extra midspan deflection is compared to
   the EN amplified-imperfection form

       delta_add = delta0 * alpha / (1 - alpha),   alpha = P / P_cr,

   which for alpha = 0.8 amplifies the seeded bow four-fold.

2. Spring-mass cantilever. A massless steel cantilever (a material with no
   density, so the beam carries no distributed mass) with a lumped point mass at
   its free tip is the textbook single-degree-of-freedom oscillator: the tip
   bending stiffness is k = 3EI/L^3, so

       f = sqrt(3EI / (m L^3)) / (2 pi).

   Modal analysis of the point-mass model is compared to that closed form.

Run after building the extension:
    cmake --build build --target _keta   (then PYTHONPATH=build/python)
    python examples/python/imperfect_column.py
"""

import math

import keta as k

# ------------------------------------------------------------------ imperfect
L = 4.0            # column length (m)
N_ELEMS = 10       # mesh density (buckling + P-Delta converge with a few)
B, H = 0.06, 0.12  # rectangular section: b along local y, h along local z (m)
ALPHA = 0.8        # working load as a fraction of the Euler critical load


def build_column(axial):
    """A pinned-pinned steel column along +X, meshed into N_ELEMS beams, under an
    axial compression `axial` (N) at the top directed toward the base. Returns
    (model, nodes)."""
    m = k.Model()
    steel = m.steel()
    rect = m.add_section_rect(B, H)
    nodes = [m.add_node(L * i / N_ELEMS, 0.0, 0.0) for i in range(N_ELEMS + 1)]
    for a, b in zip(nodes[:-1], nodes[1:]):
        m.add_beam(a, b, section=rect, material=steel)
    # Pinned-pinned: the base is a full pin that also holds twist about the axis
    # (rx) so the column has no rigid-body torsional mode; the top is a lateral
    # pin that still slides axially so the load actually compresses the member.
    m.add_support(nodes[0], dx=True, dy=True, dz=True, rx=True)
    m.add_support(nodes[-1], dy=True, dz=True)
    m.add_point_load(nodes[-1], fx=-axial, case=1)
    return m, nodes


# 1. Reference column under a unit compression: linear buckling returns lambda_1,
#    which (the reference load being 1 N) IS the Euler critical force P_cr.
m_ref, _ = build_column(1.0)
assert not m_ref.validate(), m_ref.validate()
Pcr = m_ref.buckling(case=1, n_modes=1).factors[0]
P = ALPHA * Pcr
delta0 = L / 1000.0   # the seeded erection-tolerance bow

# 2. Design column at the working compression P = 0.8 * P_cr; buckle it to recover
#    the same fundamental mode (lambda_1 should be ~ 1 / alpha = 1.25).
md, nodes = build_column(P)
buck = md.buckling(case=1, n_modes=4)
lam1 = buck.factors[0]
mid = nodes[N_ELEMS // 2]
print(buck)

# 3. Seed the fundamental buckling mode (index 0 = "mode 1") as the imperfection.
#    The mode is unit-crest (its largest free component is +1), so rescale by the
#    midspan global-Y crest to make the midspan lateral offset exactly delta0.
dy_mid = buck.mode_displacement(0, mid)[1]
scale = delta0 / dy_mid
imperfect = md.apply_imperfection(buck, mode=0, scale=scale)
print(imperfect)

# 4. Second-order (P-Delta) solve of the bowed column at P. Displacements are
#    measured FROM the imperfect geometry, so the midspan global-Y deflection is
#    exactly the load-induced ADDITIONAL bow the EN form predicts.
so = imperfect.second_order(case=1)
bow_add_fem = abs(so.displacement(mid)[1])
bow_add_en = delta0 * ALPHA / (1.0 - ALPHA)

print("\n--- 1. Imperfect pinned-pinned column (EN 1993 amplified bow) ---")
print(f"Euler P_cr (FEM buckling) : {Pcr / 1e3:10.3f} kN")
print(f"working load P = 0.8 P_cr : {P / 1e3:10.3f} kN   (lambda_1 of design model = {lam1:.4f})")
print(f"seeded bow delta0 = L/1000: {delta0 * 1e3:10.3f} mm")
print(f"amplified bow (FEM 2nd-ord): {bow_add_fem * 1e3:10.3f} mm")
print(f"delta0*alpha/(1-alpha)    : {bow_add_en * 1e3:10.3f} mm")
rel = abs(bow_add_fem - bow_add_en) / bow_add_en
print(f"relative difference       : {rel * 100:10.3f} %")
print(f"total midspan bow         : {(delta0 + bow_add_fem) * 1e3:10.3f} mm"
      f"   (delta0/(1-alpha) = {delta0 / (1 - ALPHA) * 1e3:.3f} mm)")
assert abs(lam1 - 1.0 / ALPHA) < 0.02, f"design lambda_1 should be ~1/alpha: {lam1}"
assert rel < 0.05, f"amplified bow off by {rel * 100:.2f}%"

# ------------------------------------------------------------------ point mass
E = 210e9          # Young's modulus (Pa)
IM = 8.0e-5        # section second moment of area, both axes (m^4)
LC = 2.0           # cantilever length (m)
MASS = 500.0       # lumped tip mass (kg)

mm = k.Model()
# A massless beam: a material WITHOUT a density carries no distributed mass, so
# the only inertia is the tip point mass -> a clean single-DOF oscillator.
massless = mm.add_material(e=E, nu=0.3, rho=None)
# Raw section props (Euler-Bernoulli, no shear areas => tip stiffness is exactly
# 3EI/L^3), with iy = iz = IM.
sec = mm.add_section_props(sx=0.01, iy=IM, iz=IM, ix=1.0e-6)
n1 = mm.add_node(0.0, 0.0, 0.0)
n2 = mm.add_node(LC, 0.0, 0.0)
mm.add_beam(n1, n2, section=sec, material=massless)
mm.fix(n1)
# Constrain the tip to a single bending plane: free only in dz (lateral) and ry
# (the bending rotation), so the condensed tip stiffness is 3EI/L^3 and there is
# exactly one finite spring-mass mode.
mm.add_support(n2, dx=True, dy=True, rx=True, rz=True)
pm = mm.add_point_mass(n2, MASS)   # inertia defaults to [0, 0, 0]

assert not mm.validate(), mm.validate()

modal = mm.modal(n_modes=3)
print("\n--- 2. Spring-mass cantilever (massless beam + tip point mass) ---")
print(modal)
f_fem = modal.frequencies[0]
f_hand = math.sqrt(3.0 * E * IM / (MASS * LC**3)) / (2.0 * math.pi)
rel_f = abs(f_fem - f_hand) / f_hand
print(f"tip mass                  : {MASS:10.1f} kg   (point-mass id {pm})")
print(f"f_1 (FEM modal)           : {f_fem:10.5f} Hz")
print(f"sqrt(3EI/mL^3)/2pi        : {f_hand:10.5f} Hz")
print(f"relative error            : {rel_f * 100:10.5f} %")
assert rel_f < 1e-3, f"spring-mass frequency off by {rel_f * 100:.4f}%"

print("\nOK: seeded buckling-mode imperfection amplifies as delta0*alpha/(1-alpha),")
print("    and the point-mass cantilever matches the spring-mass frequency.")
