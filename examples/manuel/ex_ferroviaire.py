# Actions ferroviaires : coefficients dynamiques, forces horizontales, balayage en
# vitesse sous HSLM-A1 d'une travée de 20 m et interaction voie-ouvrage.
import keta

Lphi = 20.0
print(f"coefficients dynamiques : phi2 = {keta.rail_phi2(Lphi):.4f}, phi3 = {keta.rail_phi3(Lphi):.4f}")
print(f"démarrage {keta.rail_traction(Lphi) / 1e3:.0f} kN, freinage LM71 {keta.rail_braking(Lphi) / 1e3:.0f} kN, "
      f"lacet {keta.rail_nosing() / 1e3:.0f} kN")
print(f"centrifuge : f = {keta.rail_centrifugal_factor(160, Lphi):.4f}, "
      f"Qtk = {keta.rail_centrifugal_force(160, 1500, 250e3, Lphi) / 1e3:.3f} kN par essieu de 250 kN")
q1, q2, e = keta.eccentric_rail(250e3)
print(f"excentricité des charges : {q1 / 1e3:.2f} / {q2 / 1e3:.2f} kN, e = {e:.4f} m")

# Travée isostatique de 20 m, 20 barres, masse 8000 kg/m3 x 1 m2
m = keta.Model()
nd = [m.add_node(float(i), 0, 0) for i in range(21)]
sec = m.add_section_props(sx=1.0, iy=0.04, iz=0.04, ix=0.08)
mat = m.add_material(e=2.0e11, nu=0.3, rho=8000.0)
barres = [m.add_beam(nd[i], nd[i + 1], sec, mat) for i in range(20)]
for i, x in enumerate(nd):
    m.add_support(x, dx=(i == 0), dy=True, dz=(i in (0, 20)), rx=True, rz=True)
f1 = m.modal(n_modes=1).frequencies[0]
print(f"première fréquence {f1:.3f} Hz ; limite inférieure n0 = {keta.rail_n0_lower(Lphi):.3f} Hz")
vitesses = [40.0 + 5.0 * i for i in range(9)]
bal = m.rail_speed_sweep(barres, "hslm_a1", vitesses, [(nd[10], 2)], zeta=0.02, n_modes=4)
for p in bal["points"]:
    print(f"  v = {p['speed']:5.1f} m/s : flèche min {p['u_min'][0] * 1e3:8.3f} mm, "
          f"accélération max {p['a_max'][0]:7.3f} m/s2")
print(f"statique {bal['static_min'][0] * 1e3:.3f} mm ; amplification dynamique {bal['dynamic_factor'][0]:.4f} ; "
      f"vitesse critique {bal['critical_speed'][0]:.1f} m/s")

# Interaction voie-ouvrage (UIC 774-3) : tablier de 60 m, appui fixe à gauche
tablier = dict(length=60, A=0.74, I=2.59, e_track=1.21, e_bearing=4.79, k_left=6e8, dT=35)
u = keta.track_bridge_uic774([tablier], rail_A=0.01534126, rail_I=6.07756e-5, rail_dT=50,
                             k_unloaded=1e7, r_unloaded=2e4, k_loaded=3e7, r_loaded=6e4,
                             train=dict(front_first=300, front_last=450, positions=2))
print(f"UIC 774-3, état thermique : contrainte du rail {u['thermal']['rail_min'] / 1e6:.2f} MPa, "
      f"effort sur l'appui fixe {u['thermal']['bearing'][0] / 1e3:.1f} kN")
print(f"enveloppe : contrainte du rail {u['envelope']['rail_min'] / 1e6:.2f} MPa, "
      f"effort sur l'appui fixe {u['envelope']['bearing'][0] / 1e3:.1f} kN")
