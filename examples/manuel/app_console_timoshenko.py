# Consoles en bois de 1, 2 et 5 m chargées en bout : part de la déformation d'effort tranchant.
import keta

E, G = 11.76e9, 0.784e9                  # bois : E / G = 15
b, h, kappa, P = 0.20, 0.50, 1.2, 50e3   # section, facteur de cisaillement, charge en bout
A, I = b * h, b * h**3 / 12

def fleche(L, cisaillement):
    m = keta.Model()
    bois = m.add_material(e=E, nu=0.3, g=G)               # G donné : aucun nu ne donne E / G = 15
    av = A / kappa if cisaillement else None              # None : rigide à l'effort tranchant
    sec = m.add_section_props(sx=A, iy=I, iz=I, ix=I, sy=av, sz=av)
    a, c = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
    m.add_beam(a, c, sec, bois)
    m.fix(a)
    m.add_point_load(c, fz=-P)
    return -m.analyze().case(1).displacement(c)[2]

print("   L     flexion    Timoshenko   forme fermée   écart   part de V")
for L in (1.0, 2.0, 5.0):
    d0, d = fleche(L, False), fleche(L, True)
    df, dv = P * L**3 / (3 * E * I), kappa * P * L / (G * A)
    ecart = 100 * (d - (df + dv)) / (df + dv)                    # en %
    print(f"{L:4.1f} m  {d0 * 1e3:7.4f} mm  {d * 1e3:7.4f} mm  {(df + dv) * 1e3:9.4f} mm  "
          f"{round(ecart, 4) + 0:6.4f} %  {100 * dv / (df + dv):6.2f} %")
    assert abs(d0 - df) < 1e-9 * df and abs(ecart) < 1e-7
