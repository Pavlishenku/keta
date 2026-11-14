# Combinaisons et enveloppes : poutre continue 2 x 10 m, charges permanentes et
# charges d'exploitation travée par travée ; combinaisons, enveloppe et extrême EN 1990.
import keta

m = keta.Model()
nd = [m.add_node(1.0 * i, 0, 0) for i in range(21)]
sec = m.add_section_props(sx=0.4, iy=0.02, iz=0.01, ix=0.01)
b = [m.add_beam(nd[i], nd[i + 1], sec, m.concrete() if i == 0 else 1) for i in range(20)]
m.add_support(nd[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(nd[10], dy=True, dz=True)
m.add_support(nd[20], dy=True, dz=True)
for bi in b:
    m.add_udl(bi, wz=-20e3, case=1)                       # G : 20 kN/m
for bi in b[:10]:
    m.add_udl(bi, wz=-15e3, case=2)                       # Q1 : travée 1
for bi in b[10:]:
    m.add_udl(bi, wz=-15e3, case=3)                       # Q2 : travée 2
for c, nom in [(1, "G"), (2, "Q travée 1"), (3, "Q travée 2")]:
    m.name_case(c, nom)
c1 = m.add_combination([(1.35, 1), (1.5, 2), (1.5, 3)], name="ELU deux travées")
c2 = m.add_combination([(1.35, 1), (1.5, 2)], name="ELU travée 1")
r = m.analyze(stations=11)
print(r)
for c in r.combo_ids:
    res = r.combo(c)
    print(f"{res} : M sur appui = {res.bar_forces(b[9])[-1][5] / 1e3:8.2f} kN.m, "
          f"M en x = 4 m = {res.bar_forces(b[3])[-1][5] / 1e3:8.2f} kN.m")

# Enveloppe programmée : G toujours présent (1,35 ou 1,0), chaque Q présent (1,5) ou absent
def efforts(cas, barre):
    return [s[5] for s in r.case(cas).bar_forces(barre)]
env = {}
for bi in b:
    G, Q1, Q2 = efforts(1, bi), efforts(2, bi), efforts(3, bi)
    vals = [[gG * g + a * q1 + c * q2 for g, q1, q2 in zip(G, Q1, Q2)]
            for gG in (1.0, 1.35) for a in (0.0, 1.5) for c in (0.0, 1.5)]
    env[bi] = [(min(v), max(v)) for v in zip(*vals)]
mmin = min(lo for bi in b for lo, _ in env[bi])
mmax = max(hi for bi in b for _, hi in env[bi])
print(f"enveloppe de My : min {mmin / 1e3:.2f} kN.m (travées), max {mmax / 1e3:.2f} kN.m (appui)")

# Extrême exact d'une réponse selon l'EN 1990 (effets des actions unitaires fournis)
e = [r.case(c).bar_forces(b[9])[-1][5] for c in (1, 2, 3)]
x = keta.combine_extreme(e, [dict(case=1, kind="permanent"),
                             dict(case=2, name="Q1", psi0=0.7),
                             dict(case=3, name="Q2", psi0=0.7)], state="uls", maximise=True)
print(f"EN 1990, ELU, max du moment sur appui : {x['value'] / 1e3:.2f} kN.m ; {x['label']} ; termes {[(c, round(f, 4)) for c, f in x['terms']]}")
x = keta.combine_extreme(e, [dict(case=1, kind="permanent"), dict(case=2, name="Q1"), dict(case=3, name="Q2")],
                         state="characteristic", maximise=True)
print(f"EN 1990, ELS caractéristique : {x['value'] / 1e3:.2f} kN.m ; termes {[(c, round(f, 4)) for c, f in x['terms']]}")
