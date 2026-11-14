# Time-history analyses: response spectrum of an accelerogram, oscillator on a lead-core
# bearing (nonlinear) and pendulum in large displacements.
import math
import keta

# 1) Response spectrum (pseudo-acceleration) of a damped sinusoidal signal of 2 Hz
dt = 0.005
acc = [3.0 * math.sin(2 * math.pi * 2.0 * k * dt) * math.exp(-0.3 * k * dt) for k in range(2001)]
T = [0.1, 0.25, 0.5, 1.0, 2.0]
Sa = keta.response_spectrum_of(acc, dt, T, zeta=0.05)
print("Sa [m/s2]:", ", ".join(f"T = {t} s: {s:.3f}" for t, s in zip(T, Sa)))

# 2) Mass of 200 t on a lead-rubber bearing, sinusoidal earthquake of 0.8 Hz
m = keta.Model()
ground, mass = m.add_node(0, 0, 0), m.add_node(0, 0, 0)
m.fix(ground)
m.add_support(mass, dy=True, dz=True, rx=True, ry=True, rz=True)
m.add_point_mass(mass, 2e5)
iso = m.add_elastic_link(ground, mass)
ts = [0.01 * k for k in range(801)]
ag = [3.0 * math.sin(2 * math.pi * 0.8 * t) if t < 5.0 else 0.0 for t in ts]
h = m.nonlinear_time_history([dict(link=iso, comp=0, law="lrb", k1=2e7, Qd=2e5, k2=1e6, C=1e5)],
                             dt=0.01, steps=800, ground_dir=0, ground=(ts, ag), observe=[(mass, 0)])
u = [x[0] for x in h["u"]]
print(f"isolator: max displacement {max(abs(x) for x in u) * 1e3:.2f} mm, "
      f"max force {max(abs(f[0]) for f in h['device_force']) / 1e3:.2f} kN")
balance = h["kinetic"] + h["strain"] + h["rayleigh_work"] + sum(h["device_work"])
print(f"input energy {h['input_work'] / 1e3:.3f} kJ = kinetic + elastic + dissipated {balance / 1e3:.3f} kJ")

# 3) Rigid pendulum of 2 m released horizontally: it passes to the other side at the half period
p = keta.Model()
a, b = p.add_node(0, 0, 0), p.add_node(2, 0, 0)
p.add_beam(a, b, p.add_section_props(1e-2, 1e-6, 1e-6, 2e-6), p.add_material(2e11, 0.3))
p.add_support(a, dx=True, dy=True, dz=True, rx=True, rz=True)
p.add_support(b, dy=True, rx=True, rz=True)
p.add_point_mass(b, 100.0)
p.add_point_load(b, fz=-981.0)
hist = p.large_displacement_time_history(1, dt=0.001, steps=1700, observe=[(b, 0)])
umin = min(x[0] for x in hist["u"])
print(f"pendulum: minimum horizontal displacement {umin:.4f} m (expected -2 L = -4 m)")
