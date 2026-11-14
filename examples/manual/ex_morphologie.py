# Simplified bridge box girder given by its outline: top slab of 12 m x 0.25 m, webs of
# 0.40 m, bottom slab of 0.30 m. The mid-line and the thicknesses are derived from the outline
# by nodes and branches, then the section is analysed as thin-walled.
import keta

contour = [
    (-6.0, 0.0), (6.0, 0.0), (6.0, -0.25), (3.0, -0.25),       # 0 to 3: top surface and right cantilever slab
    (3.0, -3.0), (-3.0, -3.0), (-3.0, -0.25), (-6.0, -0.25),   # 4 to 7: webs, soffit, left cantilever slab
    (-2.6, -0.25), (2.6, -0.25), (2.6, -2.70), (-2.6, -2.70),  # 8 to 11: cell void
    (-3.0, 0.0), (-2.6, 0.0), (2.6, 0.0), (3.0, 0.0),          # 12 to 15: top surface above the webs
]
nodes = [
    [0, 7],             # 0: left edge
    [12, 13, 6, 8],     # 1: junction of the top slab and the left web
    [14, 15, 9, 3],     # 2: junction of the top slab and the right web
    [1, 2],             # 3: right edge
    [5, 11],            # 4: bottom left corner
    [4, 10],            # 5: bottom right corner
]
branches = [(0, 1), (1, 2), (2, 3), (1, 4), (4, 5), (5, 2)]

g = keta.thinwall_from_nodes(contour, nodes, branches)
print("mid-line nodes:", [(round(y, 3), round(z, 3)) for y, z in g["nodes"]])
print("thicknesses:", [round(ta, 3) for _, _, ta, _ in g["walls"]])
p = keta.thinwall_section(g["nodes"], g["walls"])
print(f"A = {p['A']:.3f} m² (outline: {12 * 0.25 + 2 * 0.40 * 2.45 + 6 * 0.30:.3f} m², the junctions are "
      f"counted twice), zc = {p['zc']:.3f} m")
print(f"Iy = {p['Iy']:.2f} m4, J = {p['J']:.2f} m4, {p['cells']} cell")
# by hand: Bredt on the mid-line cell 5.6 x 2.725 m
Am, contour_ds_t = 5.6 * 2.725, 5.6 / 0.25 + 5.6 / 0.30 + 2 * 2.725 / 0.40
print(f"Bredt by hand: closed J = {4 * Am ** 2 / contour_ds_t:.2f} m4 (computed {p['J_closed']:.2f})")
