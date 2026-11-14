# Influence surfaces: grillage of 24 m x 9 m (5 main girders, cross-beams),
# deflection and reaction; full LM1, load groups, special vehicle and wheel set.
import keta

nx, ny, Lx, Ly = 13, 7, 24.0, 9.0
m = keta.Model()
concrete = m.concrete()
g = [[m.add_node(Lx * i / (nx - 1), Ly * j / (ny - 1), 0.0) for i in range(nx)] for j in range(ny)]
girder = m.add_section_props(sx=0.9, iy=0.35, iz=0.05, ix=0.02)     # longitudinal girders
xbeam = m.add_section_props(sx=0.3, iy=0.02, iz=0.01, ix=0.005)     # cross-beams / slab
for j in range(ny):
    for i in range(nx - 1):
        m.add_beam(g[j][i], g[j][i + 1], girder if j % 3 == 0 else xbeam, concrete)
for i in range(nx):
    for j in range(ny - 1):
        m.add_beam(g[j][i], g[j + 1][i], xbeam, concrete)
for j in range(ny):
    m.add_support(g[j][0], dx=(j == 0), dy=(j == 0), dz=True)
    m.add_support(g[j][-1], dz=True)
nodes = [nd for row in g for nd in row]                        # successive rows (order x then y)
centre = g[ny // 2][nx // 2]

srf = m.influence_surface(nodes, nx, ny, dict(kind="displacement", node=centre, comp=2))
print(f"grid {len(srf.x)} x {len(srf.y)}; ordinate at the centre {srf.at(12.0, 4.5) * 1e9:.4f} mm/MN")
print(f"deflection under 1 kN/m2: everywhere {srf.integral() * 1e6:.4f} mm, unfavourable zones "
      f"{srf.integral_neg() * 1e6:.4f} mm, rectangle [8;16]x[3;6] {srf.integral_box(8, 16, 3, 6) * 1e6:.4f} mm")
mx, mn = srf.point_extreme()
print(f"unit point load: min {mn * 1e9:.4f} mm/MN")
u = srf.udl_extreme(5e3)
print(f"distributed load of 5 kN/m2 on the unfavourable zones: {u['min'] * 1e3:.4f} mm")
lm = m.lm1_surface(nodes, nx, ny, dict(kind="displacement", node=centre, comp=2), 0.0, 9.0, minimise=True)
print(f"LM1: extreme deflection {lm['value'] * 1e3:.3f} mm; lanes:")
for v in lm["lanes"]:
    print(f"   lane {v['rank']}: edge y = {v['y_start']:.2f} m, tandem at x = {v['ts_x']:.2f} m")
gr1a = srf.load_group("gr1a", 0.0, 9.0, minimise=True)
gr4 = srf.load_group("gr4", 0.0, 9.0, minimise=True)
print(f"group gr1a: {gr1a['value'] * 1e3:.3f} mm; group gr4 (crowd): {gr4['value'] * 1e3:.3f} mm")
wheels = keta.annex_a_vehicle("1800/200", [(-1.0, 0.5), (1.0, 0.5)])  # 2 lines of wheels at +/- 1 m
v = srf.vehicle_extreme(wheels, y0=0.0, y1=9.0, minimise=True)
print(f"special vehicle 1800/200 ({len(wheels)} wheels): {v['value'] * 1e3:.3f} mm at x = {v['x']:.2f}, y = {v['y']:.2f}")
lm2 = keta.lm2_wheels()
print(f"LM2 ({lm2}): {srf.wheel_set_extreme(lm2, 48, 36)['min'] * 1e3:.3f} mm")

# Support reaction on an arbitrary cloud of nodes (triangulation of the positions in plan)
cloud = m.influence_surface_cloud(nodes, dict(kind="reaction", node=g[0][0], comp=2))
print(f"surface on a cloud: corner reaction for 100 kN at (2, 1) m = {cloud.wheels_value([(0, 0, 100e3)], 2.0, 1.0) / 1e3:.3f} kN")
