# Rupture accidentelle : poutre continue de 2 x 10 m (béton, 0,4 x 0,8 m) portée en son milieu par
# un étai ; l'étai est coupé à t = 0 depuis l'équilibre sous 20 kN/m. Flèche à mi-portée.
import keta

m = keta.Model()
mat = m.add_material(e=3.5e10, nu=0.2, rho=2500.0)
poutre = m.add_section_rect(0.4, 0.8)
etai_s = m.add_section_rect(0.3, 0.3)
nd = [m.add_node(i * 1.0, 0, 0) for i in range(21)]
barres = [m.add_beam(nd[i], nd[i + 1], poutre, mat) for i in range(20)]
sol = m.add_node(10.0, 0, -5.0)
etai = m.add_truss(sol, nd[10], section=etai_s, material=mat)
m.add_support(nd[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(nd[20], dy=True, dz=True, rx=True)
m.fix(sol)
for b in barres:
    m.add_udl(b, wz=-20e3, case=1)

T1 = 1 / m.modal(n_modes=1).frequencies[0]                 # période de la poutre intacte
h = m.rupture_history(dt=T1 / 200, steps=1600, rupture_step=0, bars=[etai],
                      cases=[(1, [0.0, 1e9], [1.0, 1.0])], observe=[(nd[10], 2)])
w0, w1 = h["static_intact"][(nd[10], 2)], h["static_severed"][(nd[10], 2)]
pic = min(h["u"][(nd[10], 2)])
print(f"flèche statique avec étai {w0 * 1e3:.3f} mm, sans étai {w1 * 1e3:.3f} mm")
print(f"flèche maximale après rupture {pic * 1e3:.3f} mm : facteur dynamique {(pic - w0) / (w1 - w0):.3f}")
