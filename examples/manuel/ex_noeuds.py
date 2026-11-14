# Génération des noeuds et des barres par des boucles Python.
import math
import keta

m = keta.Model()
# 1) Une file de 11 noeuds tous les 2 m le long de X (numéros 1 à 11)
file_x = [m.add_node(2.0 * i, 0.0, 0.0) for i in range(11)]
# 2) Un arc de cercle de rayon 15 m, 9 noeuds, dans le plan XZ
R, ouverture = 15.0, math.radians(80)
arc = []
for i in range(9):
    t = -ouverture / 2 + ouverture * i / 8
    arc.append(m.add_node(10.0 + R * math.sin(t), 8.0, R * math.cos(t) - R * math.cos(ouverture / 2)))
# 3) Une grille de 4 x 3 noeuds (numérotation par lignes)
grille = [[m.add_node(3.0 * i, 20.0 + 2.5 * j, 0.0) for i in range(4)] for j in range(3)]

acier, sec = m.steel(), m.add_section("HEA200")
barres_x = [m.add_beam(file_x[i], file_x[i + 1], sec, acier) for i in range(10)]
barres_arc = [m.add_beam(arc[i], arc[i + 1], sec, acier) for i in range(8)]
print("premier et dernier noeud de la file :", file_x[0], file_x[-1])
print("noeuds de l'arc :", arc)
print("ligne 2 de la grille :", grille[1])
print("barres de l'arc :", barres_arc)
print(m)
