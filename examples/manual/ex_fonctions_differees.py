# Functions of the time-dependent laws: EN 1992-1-1 creep and shrinkage, relaxation,
# time-dependent loss (Eq. 5.46), coded law, hereditary integral and Kelvin chain.
import math
import keta

fck, rh, h0 = 35e6, 70.0, 0.3
print(f"Ecm = {keta.ec2_ecm(fck) / 1e9:.3f} GPa")
for t0 in (7, 28, 90):
    print(f"phi(10000, {t0:2d}) = {keta.ec2_phi(10000, t0, fck, rh, h0):.4f}")
for t in (28, 365, 10000):
    print(f"shrinkage at {t:5d} days (drying from day 3): {keta.ec2_shrinkage(t, fck, rh, h0, 'N', ts=3) * 1e6:8.2f} µm/m")
print(f"age-adjusted modulus (E = 34 GPa, dphi = 2, chi = 0.8): {keta.aaem_modulus(34e9, 2.0) / 1e9:.4f} GPa")
print(f"nonlinear creep factor (|sigma| = 0.6 fcm): {keta.nonlinear_creep_factor(-0.6 * 43e6, 43e6):.4f}")

# Relaxation of a strand stressed to 0.75 fpk (class 2, rho1000 = 2.5 %) and time-dependent loss
fpk = 1860e6
for th in (1000, 5e5):
    r = keta.ec2_relaxation_ratio(th, 0.75 * fpk, fpk, 2, 2.5)
    print(f"EC2 relaxation at {th:8.0f} h: {100 * r:.3f} % of the initial stress")
print(f"BPEL relaxation at 3000 days under 1400 MPa: {keta.bpel_relaxation(3000, 1400e6, fpk, 2.5) / 1e6:.3f} MPa")
dp = keta.ec2_deferred_loss(-3e-4, 195e9, 60e6, 34e9, 2.0, 8e6, 3e-3, 0.5, 0.05, 0.3)
print(f"time-dependent loss (Eq. 5.46): {dp / 1e6:.3f} MPa")

# Coded law: time functions of a code law
law = keta.concrete_law("EC2", fc=fck, rh=rh, h_mm=300, cement="N")
print(f"EC2 law: E(28) = {law.E(28) / 1e9:.3f} GPa, E_ref(28) = {law.E_ref(28) / 1e9:.3f} GPa, "
      f"phi(10000, 28) = {law.phi(10000, 28):.4f}, J(10000, 28) = {law.compliance(10000, 28) * 1e12:.3f} 1/TPa")
t = keta.log_time_grid(28, 0.01, 1e4, 60)
eps = keta.volterra_strain(law, t, [-10e6] * len(t))
print(f"constant stress -10 MPa: final strain {eps[-1] * 1e6:.2f} µm/m "
      f"(= sigma J: {-10e6 * law.compliance(t[-1], 28) * 1e6:.2f})")
R = keta.relaxation_function(law, 28.0, t)
print(f"relaxation under a unit strain imposed at 28 days: R(28) = {R[0] / 1e9:.3f} GPa, "
      f"R(10028) = {R[-1] / 1e9:.3f} GPa")
ch = keta.fit_kelvin_chain(law, 28.0)
print(f"fitted Kelvin chain: {len(ch['tau'])} units, maximum deviation {100 * ch['max_rel_error']:.3f} %")
