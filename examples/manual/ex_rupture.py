# Accidental rupture: continuous beam of 2 x 10 m (concrete, 0.4 x 0.8 m) supported at its middle by
# a prop; the prop is cut at t = 0 from the equilibrium under 20 kN/m. Midspan deflection.
import keta

m = keta.Model()
mat = m.add_material(e=3.5e10, nu=0.2, rho=2500.0)
beam_sec = m.add_section_rect(0.4, 0.8)
prop_sec = m.add_section_rect(0.3, 0.3)
nd = [m.add_node(i * 1.0, 0, 0) for i in range(21)]
bars = [m.add_beam(nd[i], nd[i + 1], beam_sec, mat) for i in range(20)]
ground = m.add_node(10.0, 0, -5.0)
prop = m.add_truss(ground, nd[10], section=prop_sec, material=mat)
m.add_support(nd[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(nd[20], dy=True, dz=True, rx=True)
m.fix(ground)
for b in bars:
    m.add_udl(b, wz=-20e3, case=1)

T1 = 1 / m.modal(n_modes=1).frequencies[0]                 # period of the intact beam
h = m.rupture_history(dt=T1 / 200, steps=1600, rupture_step=0, bars=[prop],
                      cases=[(1, [0.0, 1e9], [1.0, 1.0])], observe=[(nd[10], 2)])
w0, w1 = h["static_intact"][(nd[10], 2)], h["static_severed"][(nd[10], 2)]
peak = min(h["u"][(nd[10], 2)])
print(f"static deflection with the prop {w0 * 1e3:.3f} mm, without the prop {w1 * 1e3:.3f} mm")
print(f"maximum deflection after rupture {peak * 1e3:.3f} mm: dynamic factor {(peak - w0) / (w1 - w0):.3f}")
