# Global activation: the same structure analysed with and without its intermediate support.
import keta

L, q = 6.0, 10e3
m = keta.Model()
n = [m.add_node(L * i, 0, 0) for i in range(3)]
s = m.add_section("IPE450")
b = [m.add_beam(n[i], n[i + 1], s, m.steel()) for i in range(2)]
m.add_support(n[0], dx=True, dy=True, dz=True, rx=True, rz=True)
m.add_support(n[1], dz=True)
m.add_support(n[2], dy=True, dz=True, rz=True)
for x in b:
    m.add_udl(x, wz=-q)
r = dict(m.analyze().case(1).reactions())
print(f"continuous beam: central support {r[n[1]][2] / 1e3:.2f} kN (5qL/4 = {1.25 * q * L / 1e3:.2f})")
m.set_active(supports=[n[0], n[2]])
r = dict(m.analyze().case(1).reactions())
print(f"central support deactivated: end support {r[n[0]][2] / 1e3:.2f} kN (qL = {q * L / 1e3:.2f})")
m.set_active()                                             # back to the full structure
print("reactivated:", round(dict(m.analyze().case(1).reactions())[n[1]][2] / 1e3, 2), "kN")
