# Sol par zones et sol pseudo-plastique : semelle quasi rigide de 10 m chargée de 1 MN
# avec un excentrement, sur un sol sans traction puis sur un sol à pression plafonnée.
import keta

P, L, kz = 1e6, 10.0, 1e7                      # charge [N], longueur [m], module x largeur [N/m par m]

def semelle(e, plastique):
    m = keta.Model()
    c = m.add_material(e=2e11, nu=0.3, rho=0.0)
    s = m.add_section_props(sx=1.0, iy=1e4, iz=1e4, ix=1e4)          # quasi rigide
    xs = sorted(set([L * i / 10 for i in range(11)] + [L / 2 + e]))
    n = [m.add_node(x, 0, 0) for x in xs]
    b = [m.add_beam(n[i], n[i + 1], s, c) for i in range(len(n) - 1)]
    m.add_support(n[0], dx=True, dy=True, rx=True, rz=True)
    m.add_point_load(n[xs.index(L / 2 + e)], fz=-P, case=1)
    r = m.pse({bi: dict(zones=[dict(kfz=kz)], plastic=plastique) for bi in b}, sub_per_zone=20)
    return n, b, r.case(1)

# 1) Sol sans traction (p = k w <= 0 en compression : seul le côté 'neg' subsiste), e = 3 m > L/6
n, b, cas = semelle(3.0, dict(side="neg"))
c = 3 * (L / 2 - 3.0)
print(f"sans traction : {cas.iterations} itérations, plastification = {cas.plastic}")
print(f"  tassement à l'about chargé {cas.displacement(n[-1])[2] * 1e3:.3f} mm "
      f"(contact sur c = {c:.1f} m, pression maximale 2P/c : {-2 * P / c / kz * 1e3:.3f} mm)")
print(f"  soulèvement de l'about opposé {cas.displacement(n[0])[2] * 1e3:.3f} mm")
sol = cas.soil(b[-1])
print(f"  dernière barre : résultante du sol {sol['total_global'][2] / 1e3:.2f} kN ; segments (composantes fx, mx, fy, fz : fz) :")
for s in sol["segments"][::5]:
    print(f"     [{s['x0']:.2f} ; {s['x1']:.2f}] m  état {s['state'][3]:8}  p moyen {s['p_mean'][3] / 1e3:9.2f} kN/m")

etats = sorted({s["state"][3] for bi in b for s in cas.soil(bi)["segments"]})
print("  états du sol rencontrés le long de la semelle :", etats)

# 2) Sol à pression plafonnée à 160 kN/m, e = 1,5 m
n, b, cas = semelle(1.5, dict(ffz=1.6e5))
segs = [s for bi in b for s in cas.soil(bi)["segments"]]
print(f"plafonné : résultante {sum(cas.soil(bi)['total_global'][2] for bi in b) / 1e3:.2f} kN, "
      f"{sum(1 for s in segs if s['state'][3] == 'capped')} segments au palier sur {len(segs)}")
