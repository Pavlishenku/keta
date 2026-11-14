# Programming the data: launching of a 90 m deck by incremental launching, on a casting bed
# (supports at x = -90, -60, -30 and 0 m) then on piers at x = 30, 60 and 90 m. For each launching
# position, a model is built, analysed and post-processed; the envelope is formed.
import csv
import io
import keta

L, mesh_step, q = 90.0, 2.5, -150e3
fixed_supports = [-90.0, -60.0, -30.0, 0.0, 30.0, 60.0, 90.0]

def deck(advance):
    """Model of the deck whose front end is at x = advance."""
    m = keta.Model()
    concrete = m.add_material(e=36e9, nu=0.2, rho=0.0)
    sec = m.add_section_props(sx=6.0, iy=5.0, iz=40.0, ix=8.0)
    xs = [advance - L + mesh_step * i for i in range(int(L / mesh_step) + 1)]
    nd = {round(x, 6): m.add_node(x, 0, 0) for x in xs}
    keys = sorted(nd)
    bars = [m.add_beam(nd[keys[i]], nd[keys[i + 1]], sec, concrete) for i in range(len(keys) - 1)]
    under = [x for x in fixed_supports if keys[0] <= x <= keys[-1]]
    for k, x in enumerate(under):
        m.add_support(nd[x], dx=(k == 0), dy=True, dz=True, rx=True)
    for b in bars:
        m.add_udl(b, wz=q)
    return m, bars, under

table = io.StringIO()                       # CSV table (here in memory; open("launching.csv", "w") otherwise)
writer = csv.writer(table, delimiter=";")
writer.writerow(["advance [m]", "supports", "M min [kN.m]", "M max [kN.m]"])
env_min, env_max = (0.0, None), (0.0, None)
for advance in [5.0 * i for i in range(19)]:
    m, bars, under = deck(advance)
    case_ = m.analyze(stations=3).case(1)
    M = [s[5] for b in bars for s in case_.bar_forces(b)]
    writer.writerow([advance, len(under), round(min(M) / 1e3, 1) + 0.0, round(max(M) / 1e3, 1) + 0.0])
    if min(M) < env_min[0]:
        env_min = (min(M), advance)
    if max(M) > env_max[0]:
        env_max = (max(M), advance)
print(table.getvalue().strip())
print(f"envelope: span moment {env_min[0] / 1e3:.1f} kN.m (advance {env_min[1]} m), "
      f"moment over support {env_max[0] / 1e3:.1f} kN.m (advance {env_max[1]} m)")
print(f"check: cantilever of 25 m, q L²/2 = {-q * 25**2 / 2 / 1e3:.1f} kN.m")
