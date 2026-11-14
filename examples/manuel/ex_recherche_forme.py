# Recherche de forme : poutre bi-appuyée de 20 m (20 kN/m) tenue à mi-portée par un hauban à 45° ;
# longueur à vide du hauban qui annule la flèche à mi-portée (effort vertical 5 q L / 8), puis la même
# chose par le défaut de longueur d'un étai.
import math
import keta

m = keta.Model()
mat = m.add_material(e=2e11, nu=0.3, rho=0.0)
tablier = m.add_section_props(sx=0.01, iy=0.04, iz=0.04, ix=0.01)        # EA = 2e9 N
n = [m.add_node(i * 1.0, 0, 0) for i in range(21)]
for i in range(20):
    m.add_udl(m.add_beam(n[i], n[i + 1], tablier, mat), wz=-2e4, case=1)
m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(n[20], dy=True, dz=True, rx=True)
pylone = m.add_node(0, 0, 10)
m.fix(pylone)
hauban = m.add_cable(pylone, n[10], ea=1e8, w=0.0, h=1e5)                 # longueur à vide : l'inconnue

r = m.tune_cables(1, [hauban], [dict(id=n[10], comp=2, value=0.0)], tol=1e-10)
F = 5 * 2e4 * 20 / 8
H = F / (1 + F / 2e9)
print(f"L0 du hauban {r['x'][0]:.6f} m en {r['iterations']} itérations ; flèche à mi-portée "
      f"{r['displacement'][n[10]][2] * 1e3:.2e} mm")
print(f"à la main : effort vertical 5qL/8 = {F / 1e3:.1f} kN, tension {math.hypot(H, F) / 1e3:.2f} kN")
