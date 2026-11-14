# Thermal actions: uniform variation on a bar clamped at both ends, gradient on a cantilever.
import keta

E, alpha, A, h = 210e9, 1.2e-5, 0.02, 0.5
m = keta.Model()
steel = m.add_material(e=E, nu=0.3, alpha_temp=alpha)
sec = m.add_section_props(sx=A, iy=8e-4, iz=2e-4, ix=1e-4)
a, b = m.add_node(0, 0, 0), m.add_node(10, 0, 0)
bar = m.add_beam(a, b, sec, steel)
m.fix(a)
m.fix(b)
m.add_thermal(bar, dt=30.0, case=1)                        # +30 K uniform
N = m.analyze().case(1).bar_forces(bar)[5][1]
print(f"fixed-fixed, +30 K: N = {N / 1e3:.3f} kN (closed form -E A alpha dT = {-E * A * alpha * 30 / 1e3:.3f} kN)")

m = keta.Model()
steel = m.add_material(e=E, nu=0.3, alpha_temp=alpha)
sec = m.add_section_props(sx=A, iy=8e-4, iz=2e-4, ix=1e-4)
a, b = m.add_node(0, 0, 0), m.add_node(10, 0, 0)
bar = m.add_beam(a, b, sec, steel)
m.fix(a)
m.add_thermal(bar, dt=0.0, gz=20.0 / h, case=1)            # 20 K difference between fibres z = -h/2 and +h/2
ux, uy, uz, rx, ry, rz = m.analyze().case(1).displacement(b)
k = alpha * 20.0 / h                                        # free curvature
print(f"cantilever, gradient 40 K/m: uz = {uz * 1e3:.3f} mm, ry = {ry * 1e3:.4f} mrad")
print(f"closed forms: |uz| = k L²/2 = {k * 100 / 2 * 1e3:.3f} mm, |ry| = k L = {k * 10 * 1e3:.4f} mrad")
