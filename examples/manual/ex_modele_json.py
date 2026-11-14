# General commands: description, validation, saving and reading back of a model.
import json
import keta

m = keta.Model()
n = [m.add_node(5.0 * i, 0.0, 0.0) for i in range(3)]
s = m.add_section_rect(0.30, 0.60)
c = m.concrete()
b = [m.add_beam(n[i], n[i + 1], s, c) for i in range(2)]
for nd in n:
    m.add_support(nd, dy=True, dz=True, rx=True)
m.add_support(n[0], dx=True)
m.add_udl(b[0], wz=-20e3, case=1)
m.name_case(1, "span 1 loaded")

print(m)                                   # model summary
print("nodes:", m.node_count, " bars:", m.bar_count)
print("problems detected:", m.validate())

text = m.to_json()                         # full model in JSON
d = json.loads(text)
print("JSON sections:", sorted(k for k in d if d[k]))
m2 = keta.Model.from_json(text)            # reading back
r1 = m.analyze().case(1).displacement(n[1])[4]
r2 = m2.analyze().case(1).displacement(n[1])[4]
print(f"rotation over support: {r1:.6e} rad (original), {r2:.6e} rad (read back)")
assert r1 == r2

script = m.export_script()                 # the model as a Python script
print("exported script:", len(script.splitlines()), "lines; start:")
print("\n".join(script.splitlines()[:3]))
print("module version:", keta.__version__)
