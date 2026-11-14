# Analyse temporelle linéaire : échelon de force sur une console à masse de tête,
# intégration directe puis superposition modale ; combinaison de directions.
import math
import keta

E, I, L, M, F = 210e9, 2e-4, 3.0, 2000.0, 10e3
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
m.add_beam(a, b, m.add_section_props(sx=0.01, iy=I, iz=10 * I, ix=I), m.add_material(e=E, nu=0.3))
m.fix(a)
m.add_point_mass(b, M)
m.add_point_load(b, fz=F, case=1)
k = 3 * E * I / L**3
T = 2 * math.pi * math.sqrt(M / k)
dt, n = T / 200, 400                                     # deux périodes
echelon = [(1, [0.0, 10 * T], [1.0, 1.0])]               # cas 1 multiplié par 1 dès t = 0
direct = m.time_history(dt=dt, steps=n, cases=echelon, observe=[(b, 2)])
modal = m.time_history(dt=dt, steps=n, cases=echelon, modal=m.modal(n_modes=1), observe=[(b, 2)])
for nom, h in (("directe", direct), ("modale ", modal)):
    print(f"{nom} : pic = {max(h['u'][(b, 2)]) * 1e3:.4f} mm  (2F/k = {2 * F / k * 1e3:.4f} mm)")

# (l'inertie latérale décuplée fait du premier mode la flexion verticale)
# combinaison de deux directions sismiques (valeurs de pointe par composante)
x = [(1.0, 0.2, 0.0, 0.0, 0.0, 0.0)]
y = [(0.3, 0.8, 0.0, 0.0, 0.0, 0.0)]
print("SRSS :", [round(v, 4) for v in keta.combine_directional([x, y], rule="srss")[0][:2]])
print("30 % :", [round(v, 4) for v in keta.combine_directional([x, y], rule="30%")[0][:2]])
