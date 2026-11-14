# Stay tuning: two stays supporting a cantilever beam of 20 m; the stay
# tensions that bring the deflection of the two attachment points back to zero are sought.
import keta

m = keta.Model()
steel = m.steel()
foot, head = m.add_node(0, 0, 0), m.add_node(0, 0, 12)
p = [m.add_node(5.0 * i, 0, 0) for i in range(1, 5)]          # x = 5, 10, 15, 20 m
deck = m.add_section_props(sx=0.05, iy=0.02, iz=0.2, ix=0.01)
pylon = m.add_section_props(sx=0.2, iy=0.5, iz=0.5, ix=0.3)
m.add_beam(foot, p[0], deck, steel)
for i in range(3):
    m.add_beam(p[i], p[i + 1], deck, steel)
m.add_beam(foot, head, pylon, steel)
stay = m.add_section_props(sx=4e-3, iy=1e-9, iz=1e-9, ix=1e-9)
h1 = m.add_truss(head, p[1], stay, steel)                      # attachment at 10 m
h2 = m.add_truss(head, p[3], stay, steel)                      # attachment at 20 m
m.fix(foot)
for b in (1, 2, 3, 4):
    m.add_udl(b, wz=-80e3, case=1)                             # deck weight
m.name_case(1, "permanent")
c1 = m.add_unit_prestress_case(h1, P=1.0)                      # unit tension cases
c2 = m.add_unit_prestress_case(h2, P=1.0)

targets = [dict(kind="displacement", id=p[1], comp=2, value=0.0),
           dict(kind="displacement", id=p[3], comp=2, value=0.0)]
t = m.tune_linear(base=[1], unknowns=[c1, c2], targets=targets)
print("tuning tensions [kN]:           ", [round(x / 1e3, 3) for x in t["x"]])
print("deflections obtained [mm]:      ", [round(x * 1e3, 9) + 0.0 for x in t["achieved"]])
print("deflections without tuning [mm]:", [round(x * 1e3, 3) for x in t["base"]])

# check: combination of the permanent case and the weighted unit cases
cmb = m.add_combination([(1.0, 1), (t["x"][0], c1), (t["x"][1], c2)], name="tuned")
res = m.analyze().combo(cmb)
print(f"combination: uz(10 m) = {res.displacement(p[1])[2] * 1e3:.6f} mm, "
      f"N stay 1 = {res.bar_forces(h1)[0][1] / 1e3:.3f} kN, N stay 2 = {res.bar_forces(h2)[0][1] / 1e3:.3f} kN")
