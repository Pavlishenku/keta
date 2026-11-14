# Calcul au feu : courbes d'incendie, propriétés des matériaux chauds, température
# d'un profilé, puis poutre BA de 6 m exposée à l'ISO 834 en sous-face et sur les côtés.
import math
import keta

for courbe in ("iso834", "external", "hydrocarbon"):
    print(f"{courbe:12}: " + ", ".join(f"{t} min {keta.fire_gas_temperature(courbe, t):6.1f} °C"
                                       for t in (15, 30, 60, 120)))
print(f"acier à 500 °C : ky = {keta.steel_ky(500):.2f}, kE = {keta.steel_ke(500):.2f} ; "
      f"température critique (mu0 = 0,5) = {keta.steel_critical_temperature(0.5):.1f} °C")
print(f"béton siliceux à 500 °C : kc = {keta.concrete_kc(500):.2f}, eps_c1 = {keta.concrete_eps_c1(500):.4f}, "
      f"dilatation = {keta.concrete_thermal_strain(500):.5f}")
print(f"armature (classe N, formée à froid) à 500 °C : ks = {keta.rebar_ks(500):.3f}, kE = {keta.rebar_kE(500):.3f}")
T = keta.steel_temperature_unprotected("iso834", am_v=200.0, ksh=1.0, t_end_minutes=30.0, dt=5.0)
print(f"profilé non protégé Am/V = 200 m-1 : {T[-1]:.1f} °C après 30 min ({len(T)} pas de 5 s)")

def poutre(q):
    """Poutre 0,20 x 0,30 m sur 6 m, 3 HA16 en haut et en bas ; beta = pi/2 rend
    l'épaisseur de la section de feu (axe local y) verticale, face 'sup' en haut."""
    m = keta.Model()
    beton = m.add_material(e=30e9, nu=0.2, rho=0.0)
    s = m.add_section_rect(0.2, 0.3)
    n = [m.add_node(1.0 * i, 0, 0) for i in range(7)]
    b = [m.add_beam(n[i], n[i + 1], s, beton, beta=math.pi / 2) for i in range(6)]
    m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
    m.add_support(n[-1], dy=True, dz=True)
    for bi in b:
        m.add_udl(bi, wz=-q, case=1)
    lit = lambda face: dict(diameter=0.016, count=3, placement=face, cover=0.04)
    feu = [dict(bar=bi, fck=30e6, rebars=[lit("inf"), lit("sup")],
                section=dict(kind="rectangle", thickness=0.3, width=0.2,
                             exposures=dict(inf=0, pz=0, mz=0), mesh=0.02)) for bi in b]
    return m, n, feu

m, n, feu = poutre(10e3)
froid = m.analyze().case(1).displacement(n[3])[2]
for duree in (30, 45, 60):
    try:
        r = m.fire(feu, [dict(kind="iso834")], duration=duree, preload_case=1, local=True)
    except ValueError:
        print(f"ISO 834, {duree} min : pas d'équilibre, ruine de la poutre avant {duree} min")
        continue
    s = r.sections()[2]
    print(f"ISO 834, {duree} min : flèche {r.displacement(n[3])[2] * 1e3:8.2f} mm (à froid {froid * 1e3:.2f}), "
          f"{r.iterations} itérations, somme des réactions {sum(x[1][2] for x in r.reactions()) / 1e3:.3f} kN")
    print(f"   section : theta max {s['theta_max']:.0f} °C, moyenne {s['theta_mean']:.0f} °C, "
          f"armatures {s['theta_rebar_max']:.0f} °C")
    print(f"   rigidités chaudes : EA = {s['EA'] / 1e6:.1f} MN, EIy = {s['EIy'] / 1e6:.3f} MN.m2, "
          f"EIz = {s['EIz'] / 1e6:.3f} MN.m2, centre (yc, zc) = ({round(s['yc'] * 1e3, 1) + 0.0:.1f}, {round(s['zc'] * 1e3, 1) + 0.0:.1f}) mm")

# Ruine : après 1 min de feu la section est à 20 °C ; Mu = 70,8 kN.m, moment moyen des barres
# centrales 13 q / 3 -> q_u = 3 Mu / 13
q_u = 3 * 70.8e3 / 13
for f in (0.97, 1.03):
    m, n, feu = poutre(f * q_u)
    try:
        m.fire(feu, [dict(kind="iso834")], duration=1, preload_case=1, local=True)
        print(f"q = {f:.2f} q_u : équilibre trouvé")
    except ValueError:
        print(f"q = {f:.2f} q_u : pas d'équilibre (ruine)")
