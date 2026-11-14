# Précontrainte par câble : poutre isostatique de 30 m en béton (EC2), câble
# parabolique tendu à 28 jours depuis l'origine ; état final à 30 ans.
import keta

L, nb = 30.0, 12
m = keta.Model()
beton = m.add_material(e=35e9, nu=0.2, rho=2500.0)
m.set_concrete_law(beton, "EC2", fc=40e6, rh=70.0, h_mm=480.0, cement="N")
b_, h_ = 0.6, 1.6
sec = m.add_section_rect(b=b_, h=h_)
n = [m.add_node(L * i / nb, 0, 0) for i in range(nb + 1)]
barres = [m.add_beam(n[i], n[i + 1], sec, beton) for i in range(nb)]
m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(n[-1], dy=True, dz=True)
m.add_self_weight(case=1)
m.name_case(1, "poids propre")

# tracé parabolique : e_z = 0 aux ancrages, -0,60 m à mi-portée (sous le centre de gravité)
trace = [(s, 0.0, -0.60 * 4 * s * (L - s) / L**2) for s in [L * i / 40 for i in range(41)]]
Ap, fpk = 3.0e-3, 1860e6
cable = m.add_tendon(barres, P0=0.75 * fpk * Ap, profile=trace, jack="start",
                     mu=0.19, k=0.005, slip=0.006, Ep=195e9, Ap=Ap, fpk=fpk, rho1000=2.5,
                     relaxation_class=2, name="C1")

r = m.phased([
    dict(name="mise en tension", time=28, cases=[1], tension=[dict(tendon=cable)]),
    dict(name="mise en service", time=90),
    dict(name="30 ans", time=28 + 30 * 365),
], creep_along_bar=True)

A, I = b_ * h_, b_ * h_**3 / 12
mi = nb // 2
for i in range(len(r)):
    p = r.phase(i)
    tid, tendu, perte, effectif = p.tendons()[0]
    x, N, Vy, Vz, T, My, Mz = p.bar_forces(barres[mi])[0]          # section à mi-portée
    s_inf = N / A + My * (-h_ / 2) / I                              # fibre inférieure
    s_sup = N / A + My * (h_ / 2) / I                               # fibre supérieure
    print(f"{p.name:16} : P = {tendu / 1e6:.3f} MN, pertes = {perte / 1e6:.3f} MN, "
          f"effectif = {effectif / 1e6:.3f} MN")
    print(f"{'':16}   flèche = {p.displacement(n[mi])[2] * 1e3:7.2f} mm, "
          f"sigma inf = {s_inf / 1e6:6.2f} MPa, sigma sup = {s_sup / 1e6:6.2f} MPa")

# Contrôle isostatique à la mise en tension : My = -q L²/8 - P(L/2) e_z(L/2), P(L/2) = -N
p = r.phase(0)
x, N, Vy, Vz, T, My, Mz = p.bar_forces(barres[mi])[0]
q = 2500.0 * 9.80665 * A
print(f"mi-portée : P(L/2) = {-N / 1e6:.4f} MN, My = {My / 1e3:.2f} kN.m, "
      f"-qL²/8 + 0,60 P = {(-q * L**2 / 8 - N * 0.60) / 1e3:.2f} kN.m")
