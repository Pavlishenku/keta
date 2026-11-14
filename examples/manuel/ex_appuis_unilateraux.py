# Appuis unilatéraux, décollement et appui élastique à matrice pleine :
# poutre à deux travées de 5 m dont seule la première est chargée.
import math
import keta

Ls, q, EI = 5.0, 10e3, 210e9 * 1e-4

def poutre():
    m = keta.Model()
    st = m.add_material(e=210e9, nu=0.3, rho=0.0)
    s = m.add_section_props(sx=0.01, iy=1e-4, iz=1e-4, ix=1e-4)
    n = [m.add_node(Ls * i / 4, 0, 0) for i in range(9)]
    b = [m.add_beam(n[i], n[i + 1], s, st) for i in range(8)]
    m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
    m.add_support(n[4], dy=True, dz=True)
    m.add_support(n[8], dy=True, dz=True)
    for bi in b[:4]:
        m.add_udl(bi, wz=-q)
    return m, n

for option in (None, "uplift", "decol"):
    m, n = poutre()
    if option == "uplift":
        m.set_support(n[8], uplift=dict(dz=+1))      # réaction verticale >= 0 seulement
    elif option == "decol":
        m.set_support(n[8], decol=("z", +1))         # tout l'appui est libéré s'il se soulève
    cas = m.analyze().case(1)
    R = dict(cas.reactions())
    Rb = R.get(n[8], (0.0,) * 6)[2]
    print(f"{str(option):6} : réactions A = {R[n[0]][2] / 1e3:7.3f}, B = {R[n[4]][2] / 1e3:7.3f}, "
          f"C = {round(Rb / 1e3, 3) + 0.0:7.3f} kN ; uz(C) = {cas.displacement(n[8])[2] * 1e3:7.3f} mm")
print(f"soulèvement libre attendu : theta_B L = q L^3 / (24 EI) x L = {q * Ls**4 / (24 * EI) * 1e3:.3f} mm")

# Appui élastique défini par sa matrice (6 termes diagonaux), tourné de 90° autour de Z :
# le ressort « x » de l'appui agit selon Y global
m = keta.Model()
st = m.add_material(e=210e9, nu=0.3, rho=0.0)
a, b = m.add_node(0, 0, 0), m.add_node(2, 0, 0)
m.add_truss(a, b, m.add_section_props(0.01, 1e-4, 1e-4, 1e-4), st)
m.fix(a)
m.add_support(b, dz=True, rx=True, ry=True, rz=True)
m.set_support(b, stiffness=[5e6, 0, 0, 0, 0, 0], rot=math.pi / 2)
m.add_point_load(b, fx=2e3, fy=1e3)
d = m.analyze().case(1).displacement(b)
print(f"appui tourné : uy = {d[1] * 1e3:.4f} mm (1 kN / 5 MN/m = 0.2000 mm), ux = {d[0] * 1e6:.4f} µm (allongement de la barre)")
