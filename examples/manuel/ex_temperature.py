# Actions thermiques : variation uniforme sur une barre bi-encastrée, gradient sur une console.
import keta

E, alpha, A, h = 210e9, 1.2e-5, 0.02, 0.5
m = keta.Model()
acier = m.add_material(e=E, nu=0.3, alpha_temp=alpha)
sec = m.add_section_props(sx=A, iy=8e-4, iz=2e-4, ix=1e-4)
a, b = m.add_node(0, 0, 0), m.add_node(10, 0, 0)
barre = m.add_beam(a, b, sec, acier)
m.fix(a)
m.fix(b)
m.add_thermal(barre, dt=30.0, case=1)                      # +30 K uniforme
N = m.analyze().case(1).bar_forces(barre)[5][1]
print(f"bi-encastrée, +30 K : N = {N / 1e3:.3f} kN (forme fermée -E A alpha dT = {-E * A * alpha * 30 / 1e3:.3f} kN)")

m = keta.Model()
acier = m.add_material(e=E, nu=0.3, alpha_temp=alpha)
sec = m.add_section_props(sx=A, iy=8e-4, iz=2e-4, ix=1e-4)
a, b = m.add_node(0, 0, 0), m.add_node(10, 0, 0)
barre = m.add_beam(a, b, sec, acier)
m.fix(a)
m.add_thermal(barre, dt=0.0, gz=20.0 / h, case=1)          # 20 K d'écart entre fibres z = -h/2 et +h/2
ux, uy, uz, rx, ry, rz = m.analyze().case(1).displacement(b)
k = alpha * 20.0 / h                                        # courbure libre
print(f"console, gradient 40 K/m : uz = {uz * 1e3:.3f} mm, ry = {ry * 1e3:.4f} mrad")
print(f"formes fermées : |uz| = k L²/2 = {k * 100 / 2 * 1e3:.3f} mm, |ry| = k L = {k * 10 * 1e3:.4f} mrad")
