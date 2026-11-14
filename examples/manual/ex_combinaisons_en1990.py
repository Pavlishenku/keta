# EN 1990 combinations generated in the model, envelopes and concomitant values.
import keta

L = 10.0
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
p = m.add_beam(a, b, m.add_section("HEB300"), m.steel())
m.add_support(a, dx=True, dy=True, dz=True, rx=True, rz=True)
m.add_support(b, dy=True, dz=True, rz=True)
m.add_udl(p, wz=-10e3, case=1); m.name_case(1, "G")
m.add_udl(p, wz=-5e3, case=2); m.name_case(2, "traffic")
m.add_udl(p, wz=2e3, case=3); m.name_case(3, "wind")
uls = m.generate_combinations([1], [(2, "traffic"), (3, "wind")], state="uls")
qp = m.generate_combinations([1], [(2, "traffic"), (3, "wind")], state="quasi_permanent")
r = m.analyze()
for i in uls + qp:
    print(f"combination {i}: reaction at A = {dict(r.combo(i).reactions())[a][2] / 1e3:8.2f} kN")
low, high = r.envelope(combos=uls)
print(f"ULS envelope of My at midspan: [{low.bar_forces(p)[5][5] / 1e3:.2f}; {high.bar_forces(p)[5][5] / 1e3:.2f}] kN.m")
cc = r.concomitant("reaction", node=a, comp=2, combos=uls)
print(f"maximum reaction {cc['max'] / 1e3:.2f} kN ({cc['source_max']}), "
      f"concomitant My = {cc['state_max'].bar_forces(p)[5][5] / 1e3:.2f} kN.m")
