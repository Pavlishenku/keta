# Poutre-caisson de hauteur variable près d'un appui : membrure supérieure horizontale (0,60 m²),
# membrure inférieure (0,35 m²) inclinée de 8 % (la hauteur croît vers l'appui), 2,60 m entre les
# centres de gravité des membrures au droit de la coupe. Moment négatif de 30 MN.m (fibre
# supérieure tendue), effort tranchant de 4 MN, effort normal de précontrainte de -12 MN.
import keta

A, Iy, Iz = 1.45, 1.92, 3.4            # section entière au droit de la coupe
z_sup, z_inf = 1.05, -1.55             # membrures, depuis le centre de gravité
membrures = [(0.60, 0.0, z_sup, 0.0, 0.0), (0.35, 0.0, z_inf, 0.0, 0.08)]
N, Vz, My = -12e6, -4e6, 30e6

r = keta.resal_web_shear(A, Iy, Iz, membrures, N=N, Vz=Vz, My=My)
print(f"effort tranchant total {Vz / 1e6:.2f} MN, repris par les membrures {r['relief_z'] / 1e6:.3f} MN, "
      f"laissé à l'âme {r['Vz'] / 1e6:.3f} MN")
# à la main : effort de la membrure inférieure, puis sa composante verticale
N_inf = 0.35 * (N / A + My * z_inf / Iy)
print(f"à la main : N_inf = {N_inf / 1e6:.3f} MN, N_inf x 0,08 = {N_inf * 0.08 / 1e6:.3f} MN")
kN, kMy, kMz = r["factors_z"]
print(f"facteurs de la section : {kN:.5f} par N, {kMy:.5f} 1/m par My")
