# EN 1993 values: plate class, effective width, buckling reduction factor.
import math
import keta

fy = 355e6
eps = keta.ec3_epsilon(fy)
print(f"epsilon = {eps:.4f} (sqrt(235/355) = {math.sqrt(235 / 355):.4f})")
# web of a welded plate girder 1200 x 12 in pure compression, flange 400 x 20
print("web c/t = 95: class", keta.ec3_class(95.0, fy), "| flange c/t = 9.7: class",
      keta.ec3_class(9.7, fy, part="outstand"))
d = keta.ec3_plate_reduction(95.0, fy)
print(f"class 4 web: k_sigma = {d['k_sigma']:.1f}, lambda_p = {d['lambda_p']:.3f}, rho = {d['rho']:.3f}")
# HEB 300 column (h/b = 1, tf = 19 mm): curve about z, then chi for lambda = 1.1
curve = keta.ec3_buckling_curve(0.30, 0.30, 0.019, about_y=False)
print(f"buckling curve about z: {curve}, chi(1.1) = {keta.ec3_chi(1.1, curve):.4f}")
