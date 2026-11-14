# Prestressing by tendon: simply supported concrete beam of 30 m (EC2), parabolic
# tendon stressed at 28 days from the start end; final state at 30 years.
import keta

L, nb = 30.0, 12
m = keta.Model()
concrete = m.add_material(e=35e9, nu=0.2, rho=2500.0)
m.set_concrete_law(concrete, "EC2", fc=40e6, rh=70.0, h_mm=480.0, cement="N")
b_, h_ = 0.6, 1.6
sec = m.add_section_rect(b=b_, h=h_)
n = [m.add_node(L * i / nb, 0, 0) for i in range(nb + 1)]
bars = [m.add_beam(n[i], n[i + 1], sec, concrete) for i in range(nb)]
m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(n[-1], dy=True, dz=True)
m.add_self_weight(case=1)
m.name_case(1, "self-weight")

# parabolic profile: e_z = 0 at the anchorages, -0.60 m at midspan (below the centroid)
profile = [(s, 0.0, -0.60 * 4 * s * (L - s) / L**2) for s in [L * i / 40 for i in range(41)]]
Ap, fpk = 3.0e-3, 1860e6
tendon = m.add_tendon(bars, P0=0.75 * fpk * Ap, profile=profile, jack="start",
                      mu=0.19, k=0.005, slip=0.006, Ep=195e9, Ap=Ap, fpk=fpk, rho1000=2.5,
                      relaxation_class=2, name="C1")

r = m.phased([
    dict(name="stressing", time=28, cases=[1], tension=[dict(tendon=tendon)]),
    dict(name="entry in service", time=90),
    dict(name="30 years", time=28 + 30 * 365),
], creep_along_bar=True)

A, I = b_ * h_, b_ * h_**3 / 12
mid = nb // 2
for i in range(len(r)):
    p = r.phase(i)
    tid, jacked, loss, effective = p.tendons()[0]
    x, N, Vy, Vz, T, My, Mz = p.bar_forces(bars[mid])[0]           # midspan section
    s_bot = N / A + My * (-h_ / 2) / I                              # bottom fibre
    s_top = N / A + My * (h_ / 2) / I                               # top fibre
    print(f"{p.name:16}: P = {jacked / 1e6:.3f} MN, losses = {loss / 1e6:.3f} MN, "
          f"effective = {effective / 1e6:.3f} MN")
    print(f"{'':16}  deflection = {p.displacement(n[mid])[2] * 1e3:7.2f} mm, "
          f"sigma bot = {s_bot / 1e6:6.2f} MPa, sigma top = {s_top / 1e6:6.2f} MPa")

# Statically determinate check at stressing: My = -q L²/8 - P(L/2) e_z(L/2), P(L/2) = -N
p = r.phase(0)
x, N, Vy, Vz, T, My, Mz = p.bar_forces(bars[mid])[0]
q = 2500.0 * 9.80665 * A
print(f"midspan: P(L/2) = {-N / 1e6:.4f} MN, My = {My / 1e3:.2f} kN.m, "
      f"-qL²/8 + 0.60 P = {(-q * L**2 / 8 - N * 0.60) / 1e3:.2f} kN.m")
