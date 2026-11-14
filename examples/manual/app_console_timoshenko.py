# Timber cantilevers of 1, 2 and 5 m loaded at the tip: share of the shear deformation.
import keta

E, G = 11.76e9, 0.784e9                  # timber: E / G = 15
b, h, kappa, P = 0.20, 0.50, 1.2, 50e3   # section, shear factor, tip load
A, I = b * h, b * h**3 / 12

def deflection(L, shear):
    m = keta.Model()
    timber = m.add_material(e=E, nu=0.3, g=G)             # G given: no nu yields E / G = 15
    av = A / kappa if shear else None                     # None: rigid in shear
    sec = m.add_section_props(sx=A, iy=I, iz=I, ix=I, sy=av, sz=av)
    a, c = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
    m.add_beam(a, c, sec, timber)
    m.fix(a)
    m.add_point_load(c, fz=-P)
    return -m.analyze().case(1).displacement(c)[2]

print("   L     bending    Timoshenko   closed form    diff.   shear part")
for L in (1.0, 2.0, 5.0):
    d0, d = deflection(L, False), deflection(L, True)
    df, dv = P * L**3 / (3 * E * I), kappa * P * L / (G * A)
    diff = 100 * (d - (df + dv)) / (df + dv)                     # in %
    print(f"{L:4.1f} m  {d0 * 1e3:7.4f} mm  {d * 1e3:7.4f} mm  {(df + dv) * 1e3:9.4f} mm  "
          f"{round(diff, 4) + 0:6.4f} %  {100 * dv / (df + dv):6.2f} %")
    assert abs(d0 - df) < 1e-9 * df and abs(diff) < 1e-7
