# Phasage sans fluage : poutre de 2 x 6 m coulée sur un étai central, puis décintrée ;
# activations, ajout et retrait d'appui, cumul des états.
import keta

L, q = 6.0, -25e3
m = keta.Model()
mat = m.add_material(e=30e9, nu=0.2, rho=0.0)        # pas de loi différée : phasage élastique
sec = m.add_section_props(sx=0.3, iy=6.4e-3, iz=2e-3, ix=2e-3)
n = [m.add_node(L * i / 2, 0, 0) for i in range(5)]
b = [m.add_beam(n[i], n[i + 1], sec, mat) for i in range(4)]
m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(n[4], dy=True, dz=True)
for bi in b:
    m.add_udl(bi, wz=q, case=1)                      # poids du béton frais

r = m.phased([
    dict(name="étaiement", time=0, activate=b, add_supports=[dict(node=n[2], dx=False, dy=True, dz=True)]),
    dict(name="coulage", time=1, cases=[1]),
    dict(name="décintrement", time=28, remove_supports=[n[2]]),
])
for i in range(len(r)):
    p = r.phase(i)
    M = p.bar_forces(b[1])[-1][5]                    # section à mi-longueur
    R = dict(p.reactions())[n[0]][2]                 # appui d'extrémité
    dM = p.increment_bar_forces(b[1])[-1][5]
    print(f"{p.name:13} : M(milieu) = {M / 1e3:8.2f} kN.m (incrément {dM / 1e3:8.2f}), "
          f"réaction en rive = {R / 1e3:7.2f} kN, flèche = {p.displacement(n[2])[2] * 1e3:7.3f} mm")
print(f"sur étai (2 travées de 6 m) : M = -q L²/8 = {-q * L * L / 8 / 1e3:.2f} kN.m, R = 3 q L/8 = {-3 * q * L / 8 / 1e3:.2f} kN")
print(f"après décintrement (12 m)   : M = q (2L)²/8 = {q * (2 * L)**2 / 8 / 1e3:.2f} kN.m, R = q L = {-q * L / 1e3:.2f} kN")
