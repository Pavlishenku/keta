# Valeurs de l'EN 1993 : classe des parois, largeur efficace, coefficient de flambement.
import math
import keta

fy = 355e6
eps = keta.ec3_epsilon(fy)
print(f"epsilon = {eps:.4f} (racine(235/355) = {math.sqrt(235 / 355):.4f})")
# âme d'un PRS 1200 x 12 en compression pure, semelle 400 x 20
print("âme c/t = 95 : classe", keta.ec3_class(95.0, fy), "| semelle c/t = 9,7 : classe",
      keta.ec3_class(9.7, fy, part="outstand"))
d = keta.ec3_plate_reduction(95.0, fy)
print(f"âme de classe 4 : k_sigma = {d['k_sigma']:.1f}, lambda_p = {d['lambda_p']:.3f}, rho = {d['rho']:.3f}")
# poteau HEB 300 (h/b = 1, tf = 19 mm) : courbe selon z, puis chi pour lambda = 1,1
courbe = keta.ec3_buckling_curve(0.30, 0.30, 0.019, about_y=False)
print(f"courbe de flambement autour de z : {courbe}, chi(1,1) = {keta.ec3_chi(1.1, courbe):.4f}")
