# Surfaces d'influence : grillage de 24 m x 9 m (5 poutres principales, entretoises),
# flèche et réaction ; LM1 complet, groupes de charges, véhicule spécial et convoi de roues.
import keta

nx, ny, Lx, Ly = 13, 7, 24.0, 9.0
m = keta.Model()
beton = m.concrete()
g = [[m.add_node(Lx * i / (nx - 1), Ly * j / (ny - 1), 0.0) for i in range(nx)] for j in range(ny)]
pp = m.add_section_props(sx=0.9, iy=0.35, iz=0.05, ix=0.02)     # poutres longitudinales
en = m.add_section_props(sx=0.3, iy=0.02, iz=0.01, ix=0.005)    # entretoises / hourdis
for j in range(ny):
    for i in range(nx - 1):
        m.add_beam(g[j][i], g[j][i + 1], pp if j % 3 == 0 else en, beton)
for i in range(nx):
    for j in range(ny - 1):
        m.add_beam(g[j][i], g[j + 1][i], en, beton)
for j in range(ny):
    m.add_support(g[j][0], dx=(j == 0), dy=(j == 0), dz=True)
    m.add_support(g[j][-1], dz=True)
noeuds = [nd for ligne in g for nd in ligne]                   # rangées successives (ordre x puis y)
centre = g[ny // 2][nx // 2]

srf = m.influence_surface(noeuds, nx, ny, dict(kind="displacement", node=centre, comp=2))
print(f"grille {len(srf.x)} x {len(srf.y)} ; ordonnée au centre {srf.at(12.0, 4.5) * 1e9:.4f} mm/MN")
print(f"flèche sous 1 kN/m2 : partout {srf.integral() * 1e6:.4f} mm, zones défavorables "
      f"{srf.integral_neg() * 1e6:.4f} mm, rectangle [8;16]x[3;6] {srf.integral_box(8, 16, 3, 6) * 1e6:.4f} mm")
mx, mn = srf.point_extreme()
print(f"charge ponctuelle unité : min {mn * 1e9:.4f} mm/MN")
u = srf.udl_extreme(5e3)
print(f"charge répartie 5 kN/m2 sur les zones défavorables : {u['min'] * 1e3:.4f} mm")
lm = m.lm1_surface(noeuds, nx, ny, dict(kind="displacement", node=centre, comp=2), 0.0, 9.0, minimise=True)
print(f"LM1 : flèche extrême {lm['value'] * 1e3:.3f} mm ; voies :")
for v in lm["lanes"]:
    print(f"   rang {v['rank']} : bord y = {v['y_start']:.2f} m, tandem en x = {v['ts_x']:.2f} m")
gr1a = srf.load_group("gr1a", 0.0, 9.0, minimise=True)
gr4 = srf.load_group("gr4", 0.0, 9.0, minimise=True)
print(f"groupe gr1a : {gr1a['value'] * 1e3:.3f} mm ; groupe gr4 (foule) : {gr4['value'] * 1e3:.3f} mm")
roues = keta.annex_a_vehicle("1800/200", [(-1.0, 0.5), (1.0, 0.5)])   # 2 files de roues à +/- 1 m
v = srf.vehicle_extreme(roues, y0=0.0, y1=9.0, minimise=True)
print(f"véhicule spécial 1800/200 ({len(roues)} roues) : {v['value'] * 1e3:.3f} mm en x = {v['x']:.2f}, y = {v['y']:.2f}")
lm2 = keta.lm2_wheels()
print(f"LM2 ({lm2}) : {srf.wheel_set_extreme(lm2, 48, 36)['min'] * 1e3:.3f} mm")

# Réaction d'appui sur un nuage de noeuds quelconque (triangulation des positions en plan)
nuage = m.influence_surface_cloud(noeuds, dict(kind="reaction", node=g[0][0], comp=2))
print(f"surface sur nuage : réaction d'angle pour 100 kN à (2, 1) m = {nuage.wheels_value([(0, 0, 100e3)], 2.0, 1.0) / 1e3:.3f} kN")
