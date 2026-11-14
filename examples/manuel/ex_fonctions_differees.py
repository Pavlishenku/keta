# Fonctions des lois différées : fluage et retrait de l'EN 1992-1-1, relaxation,
# perte différée (éq. 5.46), loi codée, intégrale héréditaire et chaîne de Kelvin.
import math
import keta

fck, rh, h0 = 35e6, 70.0, 0.3
print(f"Ecm = {keta.ec2_ecm(fck) / 1e9:.3f} GPa")
for t0 in (7, 28, 90):
    print(f"phi(10000, {t0:2d}) = {keta.ec2_phi(10000, t0, fck, rh, h0):.4f}")
for t in (28, 365, 10000):
    print(f"retrait à {t:5d} j (séchage depuis 3 j) : {keta.ec2_shrinkage(t, fck, rh, h0, 'N', ts=3) * 1e6:8.2f} µm/m")
print(f"module ajusté (E = 34 GPa, dphi = 2, chi = 0,8) : {keta.aaem_modulus(34e9, 2.0) / 1e9:.4f} GPa")
print(f"coefficient de fluage non linéaire (|sigma| = 0,6 fcm) : {keta.nonlinear_creep_factor(-0.6 * 43e6, 43e6):.4f}")

# Relaxation d'un toron tendu à 0,75 fpk (classe 2, rho1000 = 2,5 %) et perte différée
fpk = 1860e6
for th in (1000, 5e5):
    r = keta.ec2_relaxation_ratio(th, 0.75 * fpk, fpk, 2, 2.5)
    print(f"relaxation EC2 à {th:8.0f} h : {100 * r:.3f} % de la tension initiale")
print(f"relaxation BPEL à 3000 j sous 1400 MPa : {keta.bpel_relaxation(3000, 1400e6, fpk, 2.5) / 1e6:.3f} MPa")
dp = keta.ec2_deferred_loss(-3e-4, 195e9, 60e6, 34e9, 2.0, 8e6, 3e-3, 0.5, 0.05, 0.3)
print(f"perte différée (éq. 5.46) : {dp / 1e6:.3f} MPa")

# Loi codée : fonctions du temps d'une loi réglementaire
loi = keta.concrete_law("EC2", fc=fck, rh=rh, h_mm=300, cement="N")
print(f"loi EC2 : E(28) = {loi.E(28) / 1e9:.3f} GPa, E_ref(28) = {loi.E_ref(28) / 1e9:.3f} GPa, "
      f"phi(10000, 28) = {loi.phi(10000, 28):.4f}, J(10000, 28) = {loi.compliance(10000, 28) * 1e12:.3f} 1/TPa")
t = keta.log_time_grid(28, 0.01, 1e4, 60)
eps = keta.volterra_strain(loi, t, [-10e6] * len(t))
print(f"contrainte constante -10 MPa : déformation finale {eps[-1] * 1e6:.2f} µm/m "
      f"(= sigma J : {-10e6 * loi.compliance(t[-1], 28) * 1e6:.2f})")
R = keta.relaxation_function(loi, 28.0, t)
print(f"relaxation sous déformation unité imposée à 28 j : R(28) = {R[0] / 1e9:.3f} GPa, "
      f"R(10028) = {R[-1] / 1e9:.3f} GPa")
ch = keta.fit_kelvin_chain(loi, 28.0)
print(f"chaîne de Kelvin ajustée : {len(ch['tau'])} unités, écart maximal {100 * ch['max_rel_error']:.3f} %")
