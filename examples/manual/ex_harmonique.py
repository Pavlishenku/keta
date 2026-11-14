# Harmonic and random response: beam of 20 m excited at midspan by a sinusoidal
# force of amplitude 10 kN; then power spectral density (white noise).
import math
import keta

m = keta.Model()
nd = [m.add_node(i * 1.0, 0, 0) for i in range(21)]
sec = m.add_section_props(sx=1.5, iy=0.5, iz=2.0, ix=0.3)
mat = m.add_material(e=35e9, nu=0.2, rho=2500.0)
bars = [m.add_beam(nd[i], nd[i + 1], sec, mat) for i in range(20)]
for i, x in enumerate(nd):
    m.add_support(x, dx=True, dy=True, dz=(i in (0, 20)), rx=True, rz=True)
m.add_point_load(nd[10], fz=-10e3, case=1)           # force pattern (case 1 = index 0)

f1 = m.modal(n_modes=1).frequencies[0]
w_stat = m.analyze().case(1).displacement(nd[10])[2]
freqs = [f1 * (0.5 + 0.01 * i) for i in range(101)]
h = m.harmonic(case=0, freqs=freqs, zeta=0.02, n_modes=6)
fp, ap = h.peak(nd[10], 2)
print(f"f1 = {f1:.4f} Hz; static deflection {abs(w_stat) * 1e3:.4f} mm")
print(f"peak: f = {fp:.4f} Hz, amplitude {ap * 1e3:.4f} mm, amplification {ap / abs(w_stat):.3f} "
      f"(1 / 2 zeta = {1 / 0.04:.1f} for a single oscillator)")
i = freqs.index(fp)
print(f"phase at resonance: {math.degrees(h.phase(i, nd[10], 2)):.2f} degrees")
print(f"stress at the bottom fibre at midspan: {h.stress_amplitude(i, bars[10], 0, 0.0, -0.6) / 1e6:.3f} MPa")

# Random excitation F(t) = alpha(t) x pattern, spectral density of alpha = 1 /Hz (white noise 0 to 50 Hz)
fr = [0.1 * k for k in range(1, 501)]
p = m.psd(case=0, freqs=fr, zeta=0.02, n_modes=6, psd_table=[(0.0, 1.0), (50.0, 1.0)])
fpk, spk = p.peak_disp_psd(nd[10], 2)
print(f"PSD of the deflection: peak {spk * 1e6:.4f} mm2/Hz at {fpk:.2f} Hz; standard deviation {p.disp_rms(nd[10], 2) * 1e3:.4f} mm")
print(f"standard deviation of the bottom stress at midspan: {p.stress_rms(bars[10], 0, 0.0, -0.6) / 1e6:.4f} MPa")
