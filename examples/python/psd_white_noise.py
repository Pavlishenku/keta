"""keta example: random (PSD) response of an SDOF oscillator to white-noise forcing.

The same massless-bar-plus-tip-mass SDOF as harmonic_sdof.py, but now the tip
force is a random process F(t) = alpha(t) * P whose scalar multiplier alpha has a
flat (white) power spectral density S0. For any linear output the response PSD is

    S_u(f) = |H(f)|^2 * S0

where H(f) is exactly the harmonic frequency response to a unit alpha (the load
pattern P). This script drives the random solver and checks three things:

  * the PSD identity S_u = |H|^2 * S0 holds frequency by frequency, against the
    harmonic run it wraps (retrieved cheaply with psd.harmonic());
  * the resonant peak PSD equals the closed form S0 / (2*zeta*k)^2, the square of
    the harmonic resonance amplitude;
  * the RMS response sqrt(integral S_u df) is finite and positive (a truncated
    grid under-integrates the infinite band, so this is a sanity bound).

The input PSD is either psd= (a constant white level, the default 1.0) or
psd_table=[(f, S), ...] (interpolated + clamped). As with harmonic(), case= is
the 0-based dense load-case index and the damping kwargs are zeta=/zetas=/rayleigh=.

Run after building the extension:
    cmake --build build --target _keta   (then PYTHONPATH=build/python)
    python examples/python/psd_white_noise.py
"""

import math

import keta as k

E = 210e9          # Young's modulus (Pa)
A = 0.01           # axial area (m^2)
L = 2.0            # bar length (m)
MASS = 500.0       # lumped tip mass (kg)
ZETA = 0.02        # modal damping ratio (2 % of critical)
S0 = 2.0           # white input PSD level of the load multiplier (1/Hz)

stiffness = E * A / L
w_n = math.sqrt(stiffness / MASS)
f_n = w_n / (2.0 * math.pi)

m = k.Model()
mat = m.add_material(E, nu=0.3)                     # massless bar (no density)
sec = m.add_section_props(sx=A, iy=8e-5, iz=8e-5, ix=1e-6)

n1 = m.add_node(0.0, 0.0, 0.0)
n2 = m.add_node(L, 0.0, 0.0)
m.add_beam(n1, n2, section=sec, material=mat)

m.fix(n1)
m.add_support(n2, dy=True, dz=True, rx=True, ry=True, rz=True)
m.add_point_mass(n2, MASS)
m.add_point_load(n2, fx=1.0, case=1)               # UNIT load pattern P

problems = m.validate()
assert not problems, problems

# A fine grid straddling resonance (0.5 f_n .. 1.5 f_n) so the peak is resolved.
freqs = [f_n * (0.5 + i / 1000.0) for i in range(0, 1001)]
p = m.psd(case=0, freqs=freqs, zeta=ZETA, psd=S0)
print(p)

# 1. S_u(f) = |H(f)|^2 * S0 at every sweep frequency, against the wrapped FRF.
h = p.harmonic()
worst = 0.0
for fi in range(len(freqs)):
    s_u = p.disp_psd_at(fi, n2, 0)
    amp = h.amplitude_at(fi, n2, 0)
    worst = max(worst, abs(s_u - amp * amp * S0))
print(f"\nPSD identity S_u = |H|^2 * S0 : worst abs error = {worst:.2e}")
assert worst < 1e-24, f"PSD identity broken by {worst:.2e}"

# 2. Resonant peak PSD = S0 / (2*zeta*k)^2 (the square of the FRF resonance peak).
peak_f, peak_psd = p.peak_disp_psd(n2, 0)
want_peak = S0 / (2.0 * ZETA * stiffness) ** 2
print(f"\npeak PSD at          : {peak_f:9.4f} Hz  (f_n = {f_n:.4f} Hz)")
print(f"peak S_u             : {peak_psd:.6e} m^2/Hz")
print(f"S0 / (2*zeta*k)^2    : {want_peak:.6e} m^2/Hz")
rel = abs(peak_psd - want_peak) / want_peak
print(f"relative error       : {rel * 100:.3f} %")
assert rel < 2e-3, f"peak PSD off by {rel * 100:.3f}%"

# 3. RMS response is finite and positive (a truncated-band sanity bound).
rms = p.disp_rms(n2, 0)
print(f"\nRMS displacement     : {rms * 1e6:.4f} um   (sqrt(integral S_u df) over the grid)")
assert math.isfinite(rms) and rms > 0.0, "RMS must be finite and positive"
# The RMS cannot exceed the crest of the (widely sampled) resonant response.
assert rms < math.sqrt(peak_psd * (freqs[-1] - freqs[0])), "RMS bound violated"

print("\nOK: white-noise PSD response matches |H|^2 * S0, the closed-form resonant")
print("    peak S0/(2*zeta*k)^2, and a finite positive RMS.")
