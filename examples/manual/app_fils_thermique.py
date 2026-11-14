# Three wires (copper, steel, copper) carrying a rigid beam: force Q at mid-beam and
# uniform temperature rise dT of the three wires.
import keta

Ec, ac, Es, as_ = 110e9, 16.5e-6, 210e9, 12e-6    # copper, steel: moduli [Pa], expansion [1/K]
A, L, e, Q, dT = 50e-6, 0.5, 0.25, 15e3, 30.0      # area, length, spacing, force, temperature rise

m = keta.Model()
copper = m.add_material(e=Ec, nu=0.3, alpha_temp=ac)
steel = m.add_material(e=Es, nu=0.3, alpha_temp=as_)
wire = m.add_section_props(sx=A, iy=1e-12, iz=1e-12, ix=1e-12)
top = [m.add_node(x, 0, 0) for x in (-e, 0, e)]
bottom = [m.add_node(x, 0, -L) for x in (-e, 0, e)]
wires = [m.add_truss(t, b, wire, mat) for t, b, mat in zip(top, bottom, (copper, steel, copper))]
for t, b in zip(top, bottom):
    m.fix(t)
    m.add_support(b, dx=True, dy=True)
for b in (bottom[0], bottom[2]):                      # rigid beam: same vertical displacement
    m.add_constraint([(b, 2, 1.0), (bottom[1], 2, -1.0)])
m.add_point_load(bottom[1], fz=-Q)
for w in wires:
    m.add_thermal(w, dt=dT)
case = m.analyze().case(1)

# Closed form: equilibrium sigma_s + 2 sigma_c = Q / A, equal elongations of the three wires
sc = (Q / A - Es * (ac - as_) * dT) / (2 + Es / Ec)
ss = Q / A - 2 * sc
w = -(ss / Es + as_ * dT) * L
for name, v, ref in (("stress, steel  [MPa]", case.bar_forces(wires[1])[0][1] / A / 1e6, ss / 1e6),
                     ("stress, copper [MPa]", case.bar_forces(wires[0])[0][1] / A / 1e6, sc / 1e6),
                     ("beam displacement [mm]", case.displacement(bottom[1])[2] * 1e3, w * 1e3)):
    diff = 100 * (v - ref) / ref
    print(f"{name:24s} keta {v:9.4f}  closed form {ref:9.4f}  diff. {round(diff, 4) + 0:.4f} %")
    assert abs(diff) < 1e-7
