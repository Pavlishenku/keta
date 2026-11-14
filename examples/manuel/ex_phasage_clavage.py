# Phasage : deux travées préfabriquées de 20 m posées isostatiques, puis clavées
# sur l'appui central ; redistribution par fluage du moment sur appui.
import math
import keta

L, q, E = 20.0, -50e3, 35e9
m = keta.Model()
beton = m.add_material(e=E, nu=0.2, rho=0.0)
m.set_concrete_law(beton, "COURBON", e_direct=E)   # fluage non vieillissant, phi(inf) = 2
sec = m.add_section_props(sx=0.8, iy=0.25, iz=0.1, ix=0.1)   # rigide à l'effort tranchant
n = [m.add_node(L * i / 4, 0, 0) for i in range(9)]
b = [m.add_beam(n[i], n[i + 1], sec, beton) for i in range(8)]
m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(n[4], dy=True, dz=True)
m.add_support(n[8], dy=True, dz=True)
for bi in b:
    m.add_udl(bi, wz=q, case=1)                    # poids propre + superstructures

phases = [
    dict(name="pose", time=20, place_hinges=[dict(bar=b[3], end=1, dofs="ry")]),
    dict(name="chargement", time=28, cases=[1]),
    dict(name="clavage", time=60, close_hinges=[dict(bar=b[3], end=1, dofs="ry")]),
    dict(name="1 an", time=393),
    dict(name="30 ans", time=28 + 30 * 365),
]
r = m.phased(phases, creep_along_bar=True, creep_substeps_per_decade=32)
print("phases :", r.names())
for i in range(len(r)):
    p = r.phase(i)
    Ma = p.bar_forces(b[3])[-1][5] + 0.0         # moment à l'extrémité de la barre sur appui
    print(f"{p.name:11} t = {p.time:6.0f} j : M(appui) = {abs(Ma) / 1e3:8.2f} kN.m, "
          f"flèche à mi-travée = {p.displacement(n[2])[2] * 1e3:8.3f} mm")
r0 = m.phased(phases)                            # options par défaut : fluage à mi-barre
print(f"options par défaut, 30 ans : M(appui) = {r0.phase(4).bar_forces(b[3])[-1][5] / 1e3:.2f} kN.m")

# Contrôle : loi non vieillissante, limite exacte M = Mc (phi(inf) - phi(t1 - t0)) / (1 + phi(inf))
phi = lambda d: 2 * (1 - math.exp(-0.006 * d))
Mc = -q * L * L / 8                              # moment de la poutre coulée continue
M30 = Mc * (phi(30 * 365) - phi(60 - 28)) / (1 + phi(30 * 365))
print(f"poutre coulée continue : Mc = {Mc / 1e3:.2f} kN.m ; limite par fluage : {M30 / 1e3:.2f} kN.m")
