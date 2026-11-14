# Premier calcul : console en acier IPE 300 de 4 m chargée en bout.
import json
import keta

m = keta.Model()
a = m.add_node(0.0, 0.0, 0.0)            # encastrement
b = m.add_node(4.0, 0.0, 0.0)            # extrémité libre
acier = m.steel()                        # E = 210 GPa, nu = 0,3
ipe = m.add_section("IPE300")
poutre = m.add_beam(a, b, section=ipe, material=acier)
m.fix(a)
m.add_point_load(b, fz=-10e3)            # 10 kN vers le bas (cas 1)

res = m.analyze()
cas = res.case(1)
ux, uy, uz, rx, ry, rz = cas.displacement(b)
print(f"flèche en bout     uz = {uz * 1e3:.4f} mm")
print(f"rotation en bout   ry = {ry * 1e3:.4f} mrad")
noeud, (fx, fy, fz, mx, my, mz) = cas.reactions()[0]
print(f"réaction au noeud {noeud} : Fz = {fz / 1e3:.3f} kN, My = {my / 1e3:.3f} kN.m")

# Contrôle : flexion + effort tranchant (poutre de Timoshenko), avec les
# caractéristiques effectivement retenues pour la section (modèle JSON)
p = json.loads(m.to_json())["sections"][0]["props"]
E, G = 210e9, 210e9 / (2 * 1.3)
w = -10e3 * 4.0**3 / (3 * E * p["Iy"]) - 10e3 * 4.0 / (G * p["Avz"])
print(f"forme fermée       uz = {w * 1e3:.4f} mm")
assert abs(uz - w) < 1e-9 * abs(w)
