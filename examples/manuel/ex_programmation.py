# Programmation des données : lancement d'un tablier de 90 m par poussage, sur un banc
# (appuis en x = -90, -60, -30 et 0 m) puis sur des piles en x = 30, 60 et 90 m. Pour chaque position
# d'avancement, un modèle est construit, calculé et dépouillé ; on forme l'enveloppe.
import csv
import io
import keta

L, pas_maillage, q = 90.0, 2.5, -150e3
appuis_fixes = [-90.0, -60.0, -30.0, 0.0, 30.0, 60.0, 90.0]

def tablier(avance):
    """Modèle du tablier dont l'about avant est en x = avance."""
    m = keta.Model()
    beton = m.add_material(e=36e9, nu=0.2, rho=0.0)
    sec = m.add_section_props(sx=6.0, iy=5.0, iz=40.0, ix=8.0)
    xs = [avance - L + pas_maillage * i for i in range(int(L / pas_maillage) + 1)]
    nd = {round(x, 6): m.add_node(x, 0, 0) for x in xs}
    cles = sorted(nd)
    barres = [m.add_beam(nd[cles[i]], nd[cles[i + 1]], sec, beton) for i in range(len(cles) - 1)]
    sous = [x for x in appuis_fixes if cles[0] <= x <= cles[-1]]
    for k, x in enumerate(sous):
        m.add_support(nd[x], dx=(k == 0), dy=True, dz=True, rx=True)
    for b in barres:
        m.add_udl(b, wz=q)
    return m, barres, sous

table = io.StringIO()                       # tableau CSV (ici en mémoire ; open("poussage.csv", "w") sinon)
ecrit = csv.writer(table, delimiter=";")
ecrit.writerow(["avance [m]", "appuis", "M min [kN.m]", "M max [kN.m]"])
env_min, env_max = (0.0, None), (0.0, None)
for avance in [5.0 * i for i in range(19)]:
    m, barres, sous = tablier(avance)
    cas = m.analyze(stations=3).case(1)
    M = [s[5] for b in barres for s in cas.bar_forces(b)]
    ecrit.writerow([avance, len(sous), round(min(M) / 1e3, 1) + 0.0, round(max(M) / 1e3, 1) + 0.0])
    if min(M) < env_min[0]:
        env_min = (min(M), avance)
    if max(M) > env_max[0]:
        env_max = (max(M), avance)
print(table.getvalue().strip())
print(f"enveloppe : moment de travée {env_min[0] / 1e3:.1f} kN.m (avance {env_min[1]} m), "
      f"moment sur appui {env_max[0] / 1e3:.1f} kN.m (avance {env_max[1]} m)")
print(f"contrôle : console de 25 m, q L²/2 = {-q * 25**2 / 2 / 1e3:.1f} kN.m")
