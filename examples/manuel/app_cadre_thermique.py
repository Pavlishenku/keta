# Cadre en L encastré en A, sur appui glissant en C : faces extérieures à 0 °C,
# faces intérieures à +10 °C.
import keta

a, b, h = 4.0, 0.3, 0.4                   # longueur du montant et de la traverse, section b x h
E, alpha, t_ext, t_int = 30e9, 1e-5, 0.0, 10.0
I = b * h**3 / 12

m = keta.Model()
beton = m.add_material(e=E, nu=0.2, alpha_temp=alpha)
sec = m.add_section_props(sx=1e3, iy=I, iz=I, ix=I)      # axialement rigide, comme le calcul manuel
A, B, C = m.add_node(0, 0, 0), m.add_node(0, 0, a), m.add_node(a, 0, a)
montant, traverse = m.add_beam(A, B, sec, beton), m.add_beam(B, C, sec, beton)
m.fix(A)
m.add_support(C, dz=True)
for n in (B, C):
    m.add_support(n, dy=True, rx=True, rz=True)          # cadre plan XZ
# z local : -X pour le montant, +Z pour la traverse ; face intérieure du côté z < 0 des deux barres
for barre in (montant, traverse):
    m.add_thermal(barre, dt=(t_ext + t_int) / 2, gz=-(t_int - t_ext) / h)
cas = m.analyze().case(1)

X1 = -dict(cas.reactions())[C][2]                     # réaction en C, comptée vers le bas
MB = cas.bar_forces(traverse)[0][5]                   # My > 0 : fibre supérieure (froide) tendue
sigma = MB * (h / 2) / I
X1_ref = 15 * alpha * E * I / (4 * a**2) * (3 * a / h + 1)
for nom, v, ref in (("réaction X1 en C [N]", X1, X1_ref),
                    ("moment en B [N.m]", MB, a * X1_ref),
                    ("contrainte fibre froide en B [kPa]", sigma / 1e3, a * X1_ref * h / 2e3 / I)):
    ecart = 100 * (v - ref) / ref
    print(f"{nom:35s} keta {v:10.4f}  forme fermée {ref:10.4f}  écart {round(ecart, 4) + 0:.4f} %")
    assert abs(ecart) < 1e-3
