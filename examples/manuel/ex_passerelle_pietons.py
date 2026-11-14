# Passerelle sous l'action des piétons (guide Sétra/AFGC 2006, exemple A.5.1.1.1) : poutre
# mixte bi-appuyée de 38,85 m, EI = 210e9 x 0,0292, 1456 kg/m, 2,50 m de large, xi = 0,6 %,
# classe III (foule peu dense). Accélération à la résonance du premier mode vertical.
import math
import keta

L, LARGEUR, XI = 38.85, 2.5, 0.006
mu = 1456.0 + keta.setra_pedestrian_mass(0.5) * LARGEUR      # tablier + piétons de la classe III
m = keta.Model()
mat = m.add_material(e=210e9, nu=0.3, rho=mu / 100.0)
sec = m.add_section_props(sx=100.0, iy=0.0292, iz=50.0, ix=50.0)  # A fictive : inertie de rotation négligeable
nd = [m.add_node(L * i / 20, 0, 0) for i in range(21)]
for i in range(20):
    m.add_beam(nd[i], nd[i + 1], sec, mat)
m.add_support(nd[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(nd[20], dy=True, dz=True, rx=True)

modes = m.modal(n_modes=4)
f1 = modes.frequencies[0]
plage = keta.setra_frequency_range(f1, "vertical")
print(f"f1 = {f1:.3f} Hz, plage de fréquence {plage}, cas de charge {keta.setra_load_cases('III', plage)}")
s = keta.setra_surface_load("III", 1, "vertical", f1, XI, L * LARGEUR)
print(f"N = {s['pedestrians']:.1f} piétons, n_eq = {s['equivalent']:.2f}, psi = {s['psi']:.2f}, "
      f"charge {s['amplitude']:.2f} N/m2")
q = s["amplitude"] * LARGEUR                                   # charge linéique (N/m)
m.pedestrian_load(modes, 0, q, case=1)                         # signe de la déformée du mode 1
h = m.harmonic(case_id=1, freqs=[f1], zeta=XI)
acc = (2 * math.pi * f1) ** 2 * h.amplitude(nd[10], 2)[0]
print(f"accélération à mi-portée : {acc:.3f} m/s2 (forme fermée {4 / math.pi * q / mu / (2 * XI):.3f}), "
      f"plage de confort {keta.setra_comfort_range(acc)}")
