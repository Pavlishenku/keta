# Combinaisons de l'EN 1990 générées dans le modèle, enveloppes et valeurs concomitantes.
import keta

L = 10.0
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
p = m.add_beam(a, b, m.add_section("HEB300"), m.steel())
m.add_support(a, dx=True, dy=True, dz=True, rx=True, rz=True)
m.add_support(b, dy=True, dz=True, rz=True)
m.add_udl(p, wz=-10e3, case=1); m.name_case(1, "G")
m.add_udl(p, wz=-5e3, case=2); m.name_case(2, "trafic")
m.add_udl(p, wz=2e3, case=3); m.name_case(3, "vent")
elu = m.generate_combinations([1], [(2, "traffic"), (3, "wind")], state="uls")
qp = m.generate_combinations([1], [(2, "traffic"), (3, "wind")], state="quasi_permanent")
r = m.analyze()
for i in elu + qp:
    print(f"combinaison {i} : réaction en A = {dict(r.combo(i).reactions())[a][2] / 1e3:8.2f} kN")
bas, haut = r.envelope(combos=elu)
print(f"enveloppe ELU de My à mi-portée : [{bas.bar_forces(p)[5][5] / 1e3:.2f} ; {haut.bar_forces(p)[5][5] / 1e3:.2f}] kN.m")
cc = r.concomitant("reaction", node=a, comp=2, combos=elu)
print(f"réaction maximale {cc['max'] / 1e3:.2f} kN ({cc['source_max']}), "
      f"My concomitant = {cc['state_max'].bar_forces(p)[5][5] / 1e3:.2f} kN.m")
