# Réglage de tension d'un hauban de 150 m x 80 m (EA = 2e9 N, 800 N/m) : longueur à vide pour une
# tension de 2,5 MN au bas (extrémité 1) puis au haut (extrémité 2), et raideur d'Ernst de la corde.
import math
import keta

dx, dz, EA, w = 150.0, 80.0, 2e9, 800.0
for bout in (1, 2):
    l0 = keta.catenary_l0(dx, dz, ea=EA, w=w, end_tension=2.5e6, end=bout)
    c = keta.catenary(dx, dz, l0, EA, w)
    print(f"tension 2,5 MN en {bout} : L0 = {l0:.4f} m, H = {c['H'] / 1e6:.4f} MN, "
          f"T1 = {c['T1'] / 1e6:.4f} MN, T2 = {c['T2'] / 1e6:.4f} MN")
e = keta.ernst_stiffness(ea=EA, l0=l0, dx=dx, dz=dz, w=w, h=c["H"])
print(f"Ernst : k = {e['k'] / 1e6:.4f} MN/m (EA/L0 = {e['k_elastic'] / 1e6:.4f}), chi = {e['chi']:.4f}, "
      f"note {e['note']}")
