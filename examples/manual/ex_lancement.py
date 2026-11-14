# Steel deck of 60 m (0.45 m² of steel) launched on a circular alignment of 300 m radius,
# piers every 30 m. Reactions and overhang at a few positions, envelope of the moments.
import keta

R, L, n = 300.0, 60.0, 60
m = keta.Model()
steel = m.add_material(210e9, 0.3, rho=7850.0)
section = m.add_section_props(0.45, 0.60, 2.5, 0.9)          # area, Iy, Iz, J
nodes = [m.add_node(*keta.launch_path_point(L * i / n, kind="circular", radius=R)[0]) for i in range(n + 1)]
deck = [m.add_beam(a, b, section, steel) for a, b in zip(nodes, nodes[1:])]

piers = [0.0, 30.0, 60.0]
r = m.launch_on_path(deck, piers, [20.0, 45.0, 58.0, 75.0], kind="circular", radius=R)
for p in r["positions"]:
    if not p["stable"]:
        print(f"advance {p['advance']:5.1f} m: fewer than two piers under the deck")
        continue
    reactions = ", ".join(f"{v / 1e3:7.1f}" for v in p["reaction"])
    print(f"advance {p['advance']:5.1f} m: reactions [{reactions}] kN, front overhang "
          f"{p['front_overhang']:.1f} m" + (" (uplift)" if p["uplift"] else ""))
weight = 7850.0 * 0.45 * 9.80665 * L
print(f"deck weight {weight / 1e3:.1f} kN; extreme moments {min(r['m_min']) / 1e3:.0f} and "
      f"{max(r['m_max']) / 1e3:.0f} kN.m")
