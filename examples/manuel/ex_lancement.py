# Tablier métallique de 60 m (0,45 m² d'acier) poussé sur un tracé circulaire de 300 m de rayon,
# piles tous les 30 m. Réactions et porte-à-faux à quelques positions, enveloppe des moments.
import keta

R, L, n = 300.0, 60.0, 60
m = keta.Model()
acier = m.add_material(210e9, 0.3, rho=7850.0)
section = m.add_section_props(0.45, 0.60, 2.5, 0.9)          # aire, Iy, Iz, J
noeuds = [m.add_node(*keta.launch_path_point(L * i / n, kind="circular", radius=R)[0]) for i in range(n + 1)]
tablier = [m.add_beam(a, b, section, acier) for a, b in zip(noeuds, noeuds[1:])]

piles = [0.0, 30.0, 60.0]
r = m.launch_on_path(tablier, piles, [20.0, 45.0, 58.0, 75.0], kind="circular", radius=R)
for p in r["positions"]:
    if not p["stable"]:
        print(f"avance {p['advance']:5.1f} m : moins de deux piles sous le tablier")
        continue
    reactions = ", ".join(f"{v / 1e3:7.1f}" for v in p["reaction"])
    print(f"avance {p['advance']:5.1f} m : réactions [{reactions}] kN, porte-à-faux avant "
          f"{p['front_overhang']:.1f} m" + (" (soulèvement)" if p["uplift"] else ""))
poids = 7850.0 * 0.45 * 9.80665 * L
print(f"poids du tablier {poids / 1e3:.1f} kN ; moments extrêmes {min(r['m_min']) / 1e3:.0f} et "
      f"{max(r['m_max']) / 1e3:.0f} kN.m")
