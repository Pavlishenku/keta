# Fatigue : histoire du moment à mi-travée d'une poutre continue 2 x 20 m au passage
# du véhicule FLM3, comptage rainflow et dommage de Miner (courbe à pente m = 3).
import keta

m = keta.Model()
n = [m.add_node(2.0 * i, 0, 0) for i in range(21)]
b = [m.add_beam(n[i], n[i + 1], m.add_section_props(0.5, 0.08, 0.5, 0.05), m.steel()) for i in range(20)]
m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(n[10], dy=True, dz=True)
m.add_support(n[20], dy=True, dz=True)
lig = m.influence_line(b, dict(kind="bar_force", bar=b[4], x=2.0, comp=4))   # My à x = 10 m

flm3 = keta.flm3_vehicle()
print("FLM3 :", [(d, round(p / 1e3)) for d, p in flm3["axles"]], "(décalage m, kN)")
r, M = keta.convoy_history(lig, flm3)
print(f"histoire : {len(M)} positions, M de {min(M) / 1e3:.1f} à {max(M) / 1e3:.1f} kN.m")
tp = keta.turning_points(M)
print("points de rebroussement [kN.m] :", [round(x / 1e3, 1) for x in tp])
cycles = keta.rainflow(M)
for etendue, moyenne, nombre in cycles:
    print(f"   cycle : étendue {etendue / 1e3:8.2f} kN.m, moyenne {moyenne / 1e3:8.2f} kN.m, nombre {nombre}")

W = 0.08 / 0.5                                   # module de flexion I / v [m3]
dsc, NC = 71e6, 2e6                              # catégorie de détail 71 MPa à 2e6 cycles
D = sum(nb * ((e / W) / dsc) ** 3 / NC for e, _, nb in cycles)
print(f"dommage d'un passage : {D:.3e} ; passages admissibles : {1 / D:.3e}")
