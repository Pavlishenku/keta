# Wind: EN 1991-1-4 profile, mean wind, turbulent response (spectral and time-domain),
# equivalent static load and flutter of a 60 m deck.
import math
import keta

p = keta.WindProfile(vb=26.0, z0=0.05)
for z in (10.0, 30.0, 60.0):
    print(f"z = {z:4.0f} m: vm = {p.vm(z):6.3f} m/s, Iv = {p.iv(z):.4f}, qp = {p.qp(z):7.1f} Pa")
print(f"critical vortex-shedding speed (n = 1.2 Hz, d = 2.5 m, St = 0.18): "
      f"{keta.vortex_critical_speed(1.2, 2.5):.3f} m/s")
print(f"Scruton number (m = 12 t/m, zeta = 0.5 %): {keta.scruton_number(12e3, 0.005, 1.225, 2.5):.2f}")
c = keta.theodorsen(0.1)
print(f"Theodorsen function C(0.1) = {c.real:.4f} {c.imag:+.4f} i")

# Deck 12 m wide, 60 m span, 30 m high, axis along Y, wind along X
m = keta.Model()
concrete = m.add_material(e=35e9, nu=0.2, rho=2400.0)
sec = m.add_section_props(sx=6.0, iy=2.0, iz=60.0, ix=5.0)
nd = [m.add_node(0, 5.0 * i, 30) for i in range(13)]
for i in range(12):
    m.add_beam(nd[i], nd[i + 1], sec, concrete)
m.add_support(nd[0], dx=True, dy=True, dz=True, ry=True)
m.add_support(nd[-1], dx=True, dz=True, ry=True)
mid = nd[6]
spec = dict(vb=28.0, B=12.0, zeta=0.01,
            aero=dict(alpha=[-0.5, 0.5], cd=[0.6, 0.6], cl=[-1.5, 1.7], cm=[-0.1, 0.12]))
m.add_mean_wind_case(spec, case=1)
print(f"mean wind: vertical deflection at midspan {m.analyze().case(1).displacement(mid)[2] * 1e3:.3f} mm")
rep = m.buffeting(spec, [(mid, 2)], n_modes=3)["responses"][0]
print(f"turbulence (spectral): sigma = {rep['sigma'] * 1e3:.3f} mm (quasi-static {rep['sigma_b'] * 1e3:.3f}, "
      f"resonant {rep['sigma_r'] * 1e3:.3f}), peak factor g = {rep['g']:.3f}")
sigma, g = m.lrc_peak_case(spec, (mid, 2), n_modes=3, case=2)
print(f"equivalent static load (case 2): deflection {m.analyze().case(2).displacement(mid)[2] * 1e3:.3f} mm "
      f"= g sigma = {g * sigma * 1e3:.3f} mm")
dt, res = m.temporal_buffeting(spec, [(mid, 2)], n_modes=3, n_max=12.0, n_freq=2048, seed=1)
print(f"turbulence (time domain, seed 1): step {dt:.4f} s, sigma = {res[0]['sigma'] * 1e3:.3f} mm")
# Flutter: deck flexible in torsion (J = 0.06 m4), 120 m, 24 bars
m = keta.Model()
concrete = m.add_material(e=35e9, nu=0.2, rho=2500.0)
sec = m.add_section_props(sx=6.0, iy=1.2, iz=60.0, ix=0.06)
nd = [m.add_node(0, 5.0 * i, 30) for i in range(25)]
for i in range(24):
    m.add_beam(nd[i], nd[i + 1], sec, concrete)
m.add_support(nd[0], dx=True, dy=True, dz=True, ry=True)
m.add_support(nd[-1], dx=True, dz=True, ry=True)
print("natural frequencies [Hz]:", [round(f, 4) for f in m.modal(n_modes=4).frequencies])
ok, U, w = m.flutter(half_chord=6.0, a=-0.1, n_modes=2, u_min=5, u_max=400, n_sweep=120)
print(f"flutter (Theodorsen flat plate): {'found' if ok else 'absent'}, "
      f"U = {U:.3f} m/s, f = {w / (2 * math.pi):.4f} Hz")
