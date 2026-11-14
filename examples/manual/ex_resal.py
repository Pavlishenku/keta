# Box girder of variable depth near a support: horizontal top chord (0.60 m²),
# bottom chord (0.35 m²) inclined at 8 % (the depth increases towards the support), 2.60 m between the
# centroids of the chords at the cut. Negative moment of 30 MN.m (top fibre
# in tension), shear force of 4 MN, prestressing axial force of -12 MN.
import keta

A, Iy, Iz = 1.45, 1.92, 3.4            # whole section at the cut
z_top, z_bot = 1.05, -1.55             # chords, from the centroid
chords = [(0.60, 0.0, z_top, 0.0, 0.0), (0.35, 0.0, z_bot, 0.0, 0.08)]
N, Vz, My = -12e6, -4e6, 30e6

r = keta.resal_web_shear(A, Iy, Iz, chords, N=N, Vz=Vz, My=My)
print(f"total shear force {Vz / 1e6:.2f} MN, carried by the chords {r['relief_z'] / 1e6:.3f} MN, "
      f"left to the web {r['Vz'] / 1e6:.3f} MN")
# by hand: force in the bottom chord, then its vertical component
N_bot = 0.35 * (N / A + My * z_bot / Iy)
print(f"by hand: N_bot = {N_bot / 1e6:.3f} MN, N_bot x 0.08 = {N_bot * 0.08 / 1e6:.3f} MN")
kN, kMy, kMz = r["factors_z"]
print(f"section factors: {kN:.5f} per N, {kMy:.5f} 1/m per My")
