# Réponse harmonique et aléatoire : poutre de 20 m excitée à mi-portée par une force
# sinusoïdale d'amplitude 10 kN ; puis densité spectrale de puissance (bruit blanc).
import math
import keta

m = keta.Model()
nd = [m.add_node(i * 1.0, 0, 0) for i in range(21)]
sec = m.add_section_props(sx=1.5, iy=0.5, iz=2.0, ix=0.3)
mat = m.add_material(e=35e9, nu=0.2, rho=2500.0)
barres = [m.add_beam(nd[i], nd[i + 1], sec, mat) for i in range(20)]
for i, x in enumerate(nd):
    m.add_support(x, dx=True, dy=True, dz=(i in (0, 20)), rx=True, rz=True)
m.add_point_load(nd[10], fz=-10e3, case=1)           # motif de la force (cas 1 = indice 0)

f1 = m.modal(n_modes=1).frequencies[0]
w_stat = m.analyze().case(1).displacement(nd[10])[2]
freqs = [f1 * (0.5 + 0.01 * i) for i in range(101)]
h = m.harmonic(case=0, freqs=freqs, zeta=0.02, n_modes=6)
fp, ap = h.peak(nd[10], 2)
print(f"f1 = {f1:.4f} Hz ; flèche statique {abs(w_stat) * 1e3:.4f} mm")
print(f"pic : f = {fp:.4f} Hz, amplitude {ap * 1e3:.4f} mm, amplification {ap / abs(w_stat):.3f} "
      f"(1 / 2 zeta = {1 / 0.04:.1f} pour un oscillateur simple)")
i = freqs.index(fp)
print(f"phase à la résonance : {math.degrees(h.phase(i, nd[10], 2)):.2f} degrés")
print(f"contrainte en fibre inférieure à mi-portée : {h.stress_amplitude(i, barres[10], 0, 0.0, -0.6) / 1e6:.3f} MPa")

# Excitation aléatoire F(t) = alpha(t) x motif, densité spectrale de alpha = 1 /Hz (bruit blanc 0 à 50 Hz)
fr = [0.1 * k for k in range(1, 501)]
p = m.psd(case=0, freqs=fr, zeta=0.02, n_modes=6, psd_table=[(0.0, 1.0), (50.0, 1.0)])
fpk, spk = p.peak_disp_psd(nd[10], 2)
print(f"DSP de la flèche : pic {spk * 1e6:.4f} mm2/Hz à {fpk:.2f} Hz ; écart type {p.disp_rms(nd[10], 2) * 1e3:.4f} mm")
print(f"écart type de la contrainte inférieure à mi-portée : {p.stress_rms(barres[10], 0, 0.0, -0.6) / 1e6:.4f} MPa")
