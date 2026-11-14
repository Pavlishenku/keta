"""keta example: harmonic (steady-state) forced vibration of an SDOF oscillator.

A single-degree-of-freedom oscillator is the textbook check for a frequency
response. Here it is built the honest way: a massless axial bar along +X, clamped
at one end, free only in dx at the other, with a lumped tip mass. The axial
stiffness is k = E*A/L and the only inertia is the tip mass m, so

    natural angular frequency   w_n = sqrt(k / m)
    exact FRF magnitude         |u(w)| = F / sqrt((k - m*w^2)^2 + (c*w)^2)

with the viscous coefficient c = 2*zeta*sqrt(k*m). This script drives the bar
with a harmonic tip force F over a frequency sweep and checks the response
against that closed form, three ways:

  * the whole sweep matches |u(w)| to machine precision;
  * the low-frequency limit is the static response F/k, in phase with the load;
  * at resonance the amplitude is F/(2*zeta*k) -- the static response amplified
    by the resonance factor 1/(2*zeta) -- and the phase lags by exactly 90 deg.

Damping is one modal ratio (zeta=); zetas=[...] or rayleigh=(a, b) are the other
choices. The case index is the 0-based DENSE index of the load case (case=0 is
the first case), not its user id.

Run after building the extension:
    cmake --build build --target _keta   (then PYTHONPATH=build/python)
    python examples/python/harmonic_sdof.py
"""

import math

import keta as k

E = 210e9          # Young's modulus (Pa)
A = 0.01           # axial area (m^2)
L = 2.0            # bar length (m)
MASS = 500.0       # lumped tip mass (kg)
F = 1000.0         # harmonic tip force amplitude (N)
ZETA = 0.03        # modal damping ratio (3 % of critical)

stiffness = E * A / L                     # k = E*A/L
w_n = math.sqrt(stiffness / MASS)         # natural angular frequency
f_n = w_n / (2.0 * math.pi)               # natural frequency (Hz)

m = k.Model()
# A material with NO density -> a massless bar; the only inertia is the tip mass.
mat = m.add_material(E, nu=0.3)
sec = m.add_section_props(sx=A, iy=8e-5, iz=8e-5, ix=1e-6)

n1 = m.add_node(0.0, 0.0, 0.0)
n2 = m.add_node(L, 0.0, 0.0)
m.add_beam(n1, n2, section=sec, material=mat)

# Clamp the base; leave the tip free ONLY in dx -> a pure one-DOF axial oscillator.
m.fix(n1)
m.add_support(n2, dy=True, dz=True, rx=True, ry=True, rz=True)
m.add_point_mass(n2, MASS)

# The harmonic excitation: a unit-pattern axial force at the tip, load case 1.
m.add_point_load(n2, fx=F, case=1)

problems = m.validate()
assert not problems, problems

# A sweep from 0 (static) through resonance to well beyond it. Point i = 0 is
# w = 0 (the static limit) and point i = 20 lands exactly on f_n (resonance).
freqs = [f_n * i / 20.0 for i in range(0, 41)]
h = m.harmonic(case=0, freqs=freqs, zeta=ZETA)
print(h)


def frf_exact(w):
    """Exact SDOF FRF magnitude |u(w)|."""
    c = 2.0 * ZETA * math.sqrt(stiffness * MASS)
    return F / math.sqrt((stiffness - MASS * w * w) ** 2 + (c * w) ** 2)


# 1. The whole sweep matches the closed-form FRF.
worst = 0.0
for fi, f in enumerate(freqs):
    got = h.amplitude_at(fi, n2, 0)           # dx amplitude at node 2
    want = frf_exact(2.0 * math.pi * f)
    worst = max(worst, abs(got - want) / want)
print(f"\nsweep vs exact SDOF FRF : worst relative error = {worst:.2e}")
assert worst < 1e-9, f"FRF off by {worst:.2e}"

# 2. Static limit (w -> 0): amplitude F/k, in phase with the load.
static = F / stiffness
amp0 = h.amplitude_at(0, n2, 0)
phase0 = h.phase(0, n2, 0)
print(f"\nstatic amplitude (w=0)  : {amp0 * 1e6:9.4f} um   (F/k = {static * 1e6:.4f} um)")
print(f"static phase            : {math.degrees(phase0):9.4f} deg  (expect 0)")
assert abs(amp0 - static) / static < 1e-9
assert abs(phase0) < 1e-6

# 3. Resonance: peak amplitude F/(2*zeta*k) and a 90 deg phase lag.
peak_f, peak_amp = h.peak(n2, 0)
resonant = F / (2.0 * ZETA * stiffness)
amplification = peak_amp / static
phase_res = h.phase(20, n2, 0)                # index 20 sits on f_n
print(f"\npeak at                 : {peak_f:9.4f} Hz  (f_n = {f_n:.4f} Hz)")
print(f"peak amplitude          : {peak_amp * 1e6:9.4f} um   (F/(2*zeta*k) = {resonant * 1e6:.4f} um)")
print(f"resonance amplification : {amplification:9.4f}      (1/(2*zeta) = {1.0 / (2.0 * ZETA):.4f})")
print(f"phase at resonance      : {math.degrees(phase_res):9.4f} deg  (expect -90)")
assert abs(peak_f - f_n) / f_n < 1e-9, "peak should sit on the natural frequency"
assert abs(peak_amp - resonant) / resonant < 1e-3, "peak amplitude off"
assert abs(amplification - 1.0 / (2.0 * ZETA)) / (1.0 / (2.0 * ZETA)) < 1e-3
assert abs(phase_res - (-math.pi / 2.0)) < 1e-6, "resonant phase should lag 90 deg"

print("\nOK: SDOF harmonic response matches the exact FRF, the static limit and")
print("    the 1/(2*zeta) resonance amplification.")
