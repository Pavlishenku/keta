# Vent : profil EN 1991-1-4, vent moyen, réponse turbulente (spectrale et temporelle),
# charge statique équivalente et flottement d'un tablier de 60 m.
import math
import keta

p = keta.WindProfile(vb=26.0, z0=0.05)
for z in (10.0, 30.0, 60.0):
    print(f"z = {z:4.0f} m : vm = {p.vm(z):6.3f} m/s, Iv = {p.iv(z):.4f}, qp = {p.qp(z):7.1f} Pa")
print(f"vitesse critique de détachement (n = 1,2 Hz, d = 2,5 m, St = 0,18) : "
      f"{keta.vortex_critical_speed(1.2, 2.5):.3f} m/s")
print(f"nombre de Scruton (m = 12 t/m, zeta = 0,5 %) : {keta.scruton_number(12e3, 0.005, 1.225, 2.5):.2f}")
c = keta.theodorsen(0.1)
print(f"fonction de Theodorsen C(0,1) = {c.real:.4f} {c.imag:+.4f} i")

# Tablier de 12 m de large, 60 m de portée, à 30 m de hauteur, axe selon Y, vent selon X
m = keta.Model()
beton = m.add_material(e=35e9, nu=0.2, rho=2400.0)
sec = m.add_section_props(sx=6.0, iy=2.0, iz=60.0, ix=5.0)
nd = [m.add_node(0, 5.0 * i, 30) for i in range(13)]
for i in range(12):
    m.add_beam(nd[i], nd[i + 1], sec, beton)
m.add_support(nd[0], dx=True, dy=True, dz=True, ry=True)
m.add_support(nd[-1], dx=True, dz=True, ry=True)
mi = nd[6]
spec = dict(vb=28.0, B=12.0, zeta=0.01,
            aero=dict(alpha=[-0.5, 0.5], cd=[0.6, 0.6], cl=[-1.5, 1.7], cm=[-0.1, 0.12]))
m.add_mean_wind_case(spec, case=1)
print(f"vent moyen : flèche verticale à mi-portée {m.analyze().case(1).displacement(mi)[2] * 1e3:.3f} mm")
rep = m.buffeting(spec, [(mi, 2)], n_modes=3)["responses"][0]
print(f"turbulence (spectral) : sigma = {rep['sigma'] * 1e3:.3f} mm (quasi statique {rep['sigma_b'] * 1e3:.3f}, "
      f"résonant {rep['sigma_r'] * 1e3:.3f}), facteur de pointe g = {rep['g']:.3f}")
sigma, g = m.lrc_peak_case(spec, (mi, 2), n_modes=3, case=2)
print(f"charge statique équivalente (cas 2) : flèche {m.analyze().case(2).displacement(mi)[2] * 1e3:.3f} mm "
      f"= g sigma = {g * sigma * 1e3:.3f} mm")
dt, res = m.temporal_buffeting(spec, [(mi, 2)], n_modes=3, n_max=12.0, n_freq=2048, seed=1)
print(f"turbulence (temporel, graine 1) : pas {dt:.4f} s, sigma = {res[0]['sigma'] * 1e3:.3f} mm")
# Flottement : tablier souple en torsion (J = 0,06 m4), 120 m, 24 barres
m = keta.Model()
beton = m.add_material(e=35e9, nu=0.2, rho=2500.0)
sec = m.add_section_props(sx=6.0, iy=1.2, iz=60.0, ix=0.06)
nd = [m.add_node(0, 5.0 * i, 30) for i in range(25)]
for i in range(24):
    m.add_beam(nd[i], nd[i + 1], sec, beton)
m.add_support(nd[0], dx=True, dy=True, dz=True, ry=True)
m.add_support(nd[-1], dx=True, dz=True, ry=True)
print("fréquences propres [Hz] :", [round(f, 4) for f in m.modal(n_modes=4).frequencies])
ok, U, w = m.flutter(half_chord=6.0, a=-0.1, n_modes=2, u_min=5, u_max=400, n_sweep=120)
print(f"flottement (plaque plane de Theodorsen) : {'trouvé' if ok else 'absent'}, "
      f"U = {U:.3f} m/s, f = {w / (2 * math.pi):.4f} Hz")
