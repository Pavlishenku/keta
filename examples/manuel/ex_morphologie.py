# Caisson de pont simplifié donné par son contour : hourdis supérieur de 12 m x 0,25 m, âmes de
# 0,40 m, hourdis inférieur de 0,30 m. La ligne moyenne et les épaisseurs sont tirées du contour
# par nœuds et branches, puis la section est analysée en parois minces.
import keta

contour = [
    (-6.0, 0.0), (6.0, 0.0), (6.0, -0.25), (3.0, -0.25),       # 0 à 3 : extrados et encorbellement droit
    (3.0, -3.0), (-3.0, -3.0), (-3.0, -0.25), (-6.0, -0.25),   # 4 à 7 : âmes, intrados, encorbellement gauche
    (-2.6, -0.25), (2.6, -0.25), (2.6, -2.70), (-2.6, -2.70),  # 8 à 11 : alvéole
    (-3.0, 0.0), (-2.6, 0.0), (2.6, 0.0), (3.0, 0.0),          # 12 à 15 : extrados au droit des âmes
]
noeuds = [
    [0, 7],             # 0 : rive gauche
    [12, 13, 6, 8],     # 1 : jonction du hourdis supérieur et de l'âme gauche
    [14, 15, 9, 3],     # 2 : jonction du hourdis supérieur et de l'âme droite
    [1, 2],             # 3 : rive droite
    [5, 11],            # 4 : coin bas gauche
    [4, 10],            # 5 : coin bas droit
]
branches = [(0, 1), (1, 2), (2, 3), (1, 4), (4, 5), (5, 2)]

g = keta.thinwall_from_nodes(contour, noeuds, branches)
print("nœuds de la ligne moyenne :", [(round(y, 3), round(z, 3)) for y, z in g["nodes"]])
print("épaisseurs :", [round(ta, 3) for _, _, ta, _ in g["walls"]])
p = keta.thinwall_section(g["nodes"], g["walls"])
print(f"A = {p['A']:.3f} m² (contour : {12 * 0.25 + 2 * 0.40 * 2.45 + 6 * 0.30:.3f} m², les jonctions sont "
      f"comptées deux fois), zc = {p['zc']:.3f} m")
print(f"Iy = {p['Iy']:.2f} m4, J = {p['J']:.2f} m4, {p['cells']} cellule")
# à la main : Bredt sur la cellule de ligne moyenne 5,6 x 2,725 m
Am, contour_ds_t = 5.6 * 2.725, 5.6 / 0.25 + 5.6 / 0.30 + 2 * 2.725 / 0.40
print(f"Bredt à la main : J fermé = {4 * Am ** 2 / contour_ds_t:.2f} m4 (calculé {p['J_closed']:.2f})")
