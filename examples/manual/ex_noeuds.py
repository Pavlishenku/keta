# Generation of nodes and bars by Python loops.
import math
import keta

m = keta.Model()
# 1) A line of 11 nodes every 2 m along X (numbers 1 to 11)
line_x = [m.add_node(2.0 * i, 0.0, 0.0) for i in range(11)]
# 2) A circular arc of radius 15 m, 9 nodes, in the XZ plane
R, angle = 15.0, math.radians(80)
arc = []
for i in range(9):
    t = -angle / 2 + angle * i / 8
    arc.append(m.add_node(10.0 + R * math.sin(t), 8.0, R * math.cos(t) - R * math.cos(angle / 2)))
# 3) A grid of 4 x 3 nodes (numbered row by row)
grid = [[m.add_node(3.0 * i, 20.0 + 2.5 * j, 0.0) for i in range(4)] for j in range(3)]

steel, sec = m.steel(), m.add_section("HEA200")
bars_x = [m.add_beam(line_x[i], line_x[i + 1], sec, steel) for i in range(10)]
bars_arc = [m.add_beam(arc[i], arc[i + 1], sec, steel) for i in range(8)]
print("first and last node of the line:", line_x[0], line_x[-1])
print("nodes of the arc:", arc)
print("row 2 of the grid:", grid[1])
print("bars of the arc:", bars_arc)
print(m)
