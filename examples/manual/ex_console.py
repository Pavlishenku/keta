# First analysis: 4 m IPE 300 steel cantilever loaded at its tip.
import json
import keta

m = keta.Model()
a = m.add_node(0.0, 0.0, 0.0)            # fixed support
b = m.add_node(4.0, 0.0, 0.0)            # free end
steel = m.steel()                        # E = 210 GPa, nu = 0.3
ipe = m.add_section("IPE300")
beam = m.add_beam(a, b, section=ipe, material=steel)
m.fix(a)
m.add_point_load(b, fz=-10e3)            # 10 kN downwards (case 1)

res = m.analyze()
case_ = res.case(1)
ux, uy, uz, rx, ry, rz = case_.displacement(b)
print(f"tip deflection     uz = {uz * 1e3:.4f} mm")
print(f"tip rotation       ry = {ry * 1e3:.4f} mrad")
node, (fx, fy, fz, mx, my, mz) = case_.reactions()[0]
print(f"reaction at node {node}: Fz = {fz / 1e3:.3f} kN, My = {my / 1e3:.3f} kN.m")

# Check: bending + shear (Timoshenko beam), with the properties
# actually used for the section (JSON model)
p = json.loads(m.to_json())["sections"][0]["props"]
E, G = 210e9, 210e9 / (2 * 1.3)
w = -10e3 * 4.0**3 / (3 * E * p["Iy"]) - 10e3 * 4.0 / (G * p["Avz"])
print(f"closed form        uz = {w * 1e3:.4f} mm")
assert abs(uz - w) < 1e-9 * abs(w)
