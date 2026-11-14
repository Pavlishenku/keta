# Imposed strains: uniform curvature, shortening, linear variation, dislocation.
import keta

E, A, I, L = 30e9, 0.2, 4e-3, 5.0
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
p = m.add_beam(a, b, m.add_section_props(sx=A, iy=I, iz=I, ix=I), m.add_material(e=E, nu=0.2))
m.fix(a)
k = 2e-4
m.add_prestrain(p, ky=k, case=1)                        # uniform curvature
m.add_prestrain(p, strain=-3e-4, case=2)                # uniform shortening
m.add_prestrain_linear(p, (0, 0, 0), (0, k, 0), case=3)  # curvature from 0 to k
m.add_prestrain_point(p, 0.5, ry=1e-3, case=4)          # 1 mrad kink at mid-length
r = m.analyze()
d = [r.case(c).displacement(b) for c in (1, 2, 3, 4)]
print(f"uniform curvature:  rotation {d[0][4] * 1e3:8.4f} mrad (kL = {k * L * 1e3:.4f}), "
      f"deflection {d[0][2] * 1e3:8.4f} mm (kL2/2 = {k * L * L / 2 * 1e3:.4f})")
print(f"shortening:         ux = {d[1][0] * 1e3:8.4f} mm (eps L = {-3e-4 * L * 1e3:.4f})")
print(f"linear curvature:   rotation {d[2][4] * 1e3:8.4f} mrad (kL/2 = {k * L / 2 * 1e3:.4f})")
print(f"kink at mid-length: deflection {d[3][2] * 1e3:8.4f} mm (1e-3 L/2 = {1e-3 * L / 2 * 1e3:.4f})")

# same strains on a bar clamped at both ends: forces, no displacement
m.fix(b)
st = m.analyze().case(2).bar_forces(p)[5]
print(f"restrained shortening: N = {st[1] / 1e3:.1f} kN (-E A eps = {E * A * 3e-4 / 1e3:.1f})")
