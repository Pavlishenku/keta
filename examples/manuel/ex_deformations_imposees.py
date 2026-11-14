# Déformations imposées : courbure uniforme, raccourcissement, variation linéaire, dislocation.
import keta

E, A, I, L = 30e9, 0.2, 4e-3, 5.0
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
p = m.add_beam(a, b, m.add_section_props(sx=A, iy=I, iz=I, ix=I), m.add_material(e=E, nu=0.2))
m.fix(a)
k = 2e-4
m.add_prestrain(p, ky=k, case=1)                        # courbure uniforme
m.add_prestrain(p, strain=-3e-4, case=2)                # raccourcissement uniforme
m.add_prestrain_linear(p, (0, 0, 0), (0, k, 0), case=3)  # courbure de 0 à k
m.add_prestrain_point(p, 0.5, ry=1e-3, case=4)          # pli de 1 mrad à mi-longueur
r = m.analyze()
d = [r.case(c).displacement(b) for c in (1, 2, 3, 4)]
print(f"courbure uniforme  : rotation {d[0][4] * 1e3:8.4f} mrad (kL = {k * L * 1e3:.4f}), "
      f"flèche {d[0][2] * 1e3:8.4f} mm (kL2/2 = {k * L * L / 2 * 1e3:.4f})")
print(f"raccourcissement   : ux = {d[1][0] * 1e3:8.4f} mm (eps L = {-3e-4 * L * 1e3:.4f})")
print(f"courbure linéaire  : rotation {d[2][4] * 1e3:8.4f} mrad (kL/2 = {k * L / 2 * 1e3:.4f})")
print(f"pli à mi-longueur  : flèche {d[3][2] * 1e3:8.4f} mm (1e-3 L/2 = {1e-3 * L / 2 * 1e3:.4f})")

# mêmes déformations sur une barre encastrée aux deux bouts : efforts, sans déplacement
m.fix(b)
st = m.analyze().case(2).bar_forces(p)[5]
print(f"raccourcissement empêché : N = {st[1] / 1e3:.1f} kN (-E A eps = {E * A * 3e-4 / 1e3:.1f})")
