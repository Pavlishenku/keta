# Footbridge under pedestrian action (Sétra/AFGC 2006 guide, example A.5.1.1.1): simply
# supported composite beam of 38.85 m, EI = 210e9 x 0.0292, 1456 kg/m, 2.50 m wide, xi = 0.6 %,
# class III (sparse crowd). Acceleration at resonance of the first vertical mode.
import math
import keta

L, WIDTH, XI = 38.85, 2.5, 0.006
mu = 1456.0 + keta.setra_pedestrian_mass(0.5) * WIDTH        # deck + pedestrians of class III
m = keta.Model()
mat = m.add_material(e=210e9, nu=0.3, rho=mu / 100.0)
sec = m.add_section_props(sx=100.0, iy=0.0292, iz=50.0, ix=50.0)  # fictitious A: negligible rotary inertia
nd = [m.add_node(L * i / 20, 0, 0) for i in range(21)]
for i in range(20):
    m.add_beam(nd[i], nd[i + 1], sec, mat)
m.add_support(nd[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(nd[20], dy=True, dz=True, rx=True)

modes = m.modal(n_modes=4)
f1 = modes.frequencies[0]
freq_range = keta.setra_frequency_range(f1, "vertical")
print(f"f1 = {f1:.3f} Hz, frequency range {freq_range}, load cases {keta.setra_load_cases('III', freq_range)}")
s = keta.setra_surface_load("III", 1, "vertical", f1, XI, L * WIDTH)
print(f"N = {s['pedestrians']:.1f} pedestrians, n_eq = {s['equivalent']:.2f}, psi = {s['psi']:.2f}, "
      f"load {s['amplitude']:.2f} N/m2")
q = s["amplitude"] * WIDTH                                     # distributed load (N/m)
m.pedestrian_load(modes, 0, q, case=1)                         # sign of the mode 1 shape
h = m.harmonic(case_id=1, freqs=[f1], zeta=XI)
acc = (2 * math.pi * f1) ** 2 * h.amplitude(nd[10], 2)[0]
print(f"acceleration at midspan: {acc:.3f} m/s2 (closed form {4 / math.pi * q / mu / (2 * XI):.3f}), "
      f"comfort range {keta.setra_comfort_range(acc)}")
