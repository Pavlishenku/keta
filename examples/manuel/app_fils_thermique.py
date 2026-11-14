# Trois fils (cuivre, acier, cuivre) portant une poutre rigide : force Q au milieu et
# échauffement uniforme dT des trois fils.
import keta

Ec, ac, Es, as_ = 110e9, 16.5e-6, 210e9, 12e-6    # cuivre, acier : modules [Pa], dilatation [1/K]
A, L, e, Q, dT = 50e-6, 0.5, 0.25, 15e3, 30.0      # aire, longueur, entraxe, force, échauffement

m = keta.Model()
cuivre = m.add_material(e=Ec, nu=0.3, alpha_temp=ac)
acier = m.add_material(e=Es, nu=0.3, alpha_temp=as_)
fil = m.add_section_props(sx=A, iy=1e-12, iz=1e-12, ix=1e-12)
haut = [m.add_node(x, 0, 0) for x in (-e, 0, e)]
bas = [m.add_node(x, 0, -L) for x in (-e, 0, e)]
fils = [m.add_truss(h, b, fil, mat) for h, b, mat in zip(haut, bas, (cuivre, acier, cuivre))]
for h, b in zip(haut, bas):
    m.fix(h)
    m.add_support(b, dx=True, dy=True)
for b in (bas[0], bas[2]):                            # poutre rigide : même déplacement vertical
    m.add_constraint([(b, 2, 1.0), (bas[1], 2, -1.0)])
m.add_point_load(bas[1], fz=-Q)
for f in fils:
    m.add_thermal(f, dt=dT)
cas = m.analyze().case(1)

# Forme fermée : équilibre sigma_s + 2 sigma_c = Q / A, égalité des allongements des trois fils
sc = (Q / A - Es * (ac - as_) * dT) / (2 + Es / Ec)
ss = Q / A - 2 * sc
w = -(ss / Es + as_ * dT) * L
for nom, v, ref in (("contrainte acier  [MPa]", cas.bar_forces(fils[1])[0][1] / A / 1e6, ss / 1e6),
                    ("contrainte cuivre [MPa]", cas.bar_forces(fils[0])[0][1] / A / 1e6, sc / 1e6),
                    ("déplacement de la poutre [mm]", cas.displacement(bas[1])[2] * 1e3, w * 1e3)):
    ecart = 100 * (v - ref) / ref
    print(f"{nom:30s} keta {v:9.4f}  forme fermée {ref:9.4f}  écart {round(ecart, 4) + 0:.4f} %")
    assert abs(ecart) < 1e-7
