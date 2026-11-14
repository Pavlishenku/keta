# Buckling and second order: pinned-pinned column of 6 m (HEB 200), buckling about the weak axis.
# beta = pi/2 orients the weak axis so that buckling occurs along X (see the note
# on the local axes of nearly vertical bars).
import math
import keta

L, n, P = 6.0, 12, 100e3
def column():
    m = keta.Model()
    nd = [m.add_node(0, 0, L * i / n) for i in range(n + 1)]
    sec, steel = m.add_section("HEB200"), m.steel()
    b = [m.add_beam(nd[i], nd[i + 1], sec, steel, beta=math.pi / 2) for i in range(n)]
    m.add_support(nd[0], dx=True, dy=True, dz=True, rz=True)
    m.add_support(nd[-1], dx=True, dy=True)
    m.add_point_load(nd[-1], fz=-P, case=1)
    return m, nd, b

m, nd, b = column()
bk = m.buckling(case=1, n_modes=3)
print(bk)
Iz = keta.section_props("HEB200")["Iz"]
Pcr = math.pi**2 * 210e9 * Iz / L**2
print(f"critical factors: {[round(x, 4) for x in bk.factors]}")
print(f"Euler (weak axis): Pcr / P = {Pcr / P:.4f}")
print(f"mode 1: components of the mid-height node = {[round(v, 4) + 0.0 for v in bk.mode_displacement(0, nd[n // 2])]}")

# Initial imperfection in the shape of mode 1, amplitude L/300, then second-order analysis
imperfect = m.apply_imperfection(bk, mode=0, scale=L / 300)
so = imperfect.second_order(case=1, max_iterations=50)
u = max(abs(so.displacement(x)[0]) for x in nd)
print(f"second order: {so.iterations} iterations, max additional displacement {u * 1e3:.3f} mm")
amp = 1 / (1 - P / Pcr)
print(f"theoretical amplification e0 / (Pcr/P - 1) = {L / 300 / (Pcr / P - 1) * 1e3:.3f} mm "
      f"(factor 1 / (1 - P/Pcr) = {amp:.4f})")
Mz = max(abs(s[6]) for bi in b for s in so.bar_forces(bi))
print(f"maximum second-order moment: {Mz / 1e3:.3f} kN.m (P e0 amp = {P * L / 300 * amp / 1e3:.3f} kN.m)")
