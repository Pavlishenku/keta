# Natural frequencies of a simply supported beam, bare and then carrying at mid-span a spring,
# a simple support or a point mass.
import math
import keta

E, rho, L, a = 2e11, 7800.0, 1.0, 0.05     # steel, span, square section a x a
A, I, n = a * a, a**4 / 12, 80
k_spring, mass = 1e6, 0.2                  # N/m, kg

def keta_f(attachment, material_mass=False):
    m = keta.Model()
    steel = m.add_material(e=E, nu=0.3, rho=rho if material_mass else None)
    sec = m.add_section_props(sx=A, iy=I, iz=I, ix=I)
    nd = [m.add_node(L * i / n, 0, 0) for i in range(n + 1)]
    for i in range(n):
        b = m.add_beam(nd[i], nd[i + 1], sec, steel)
        if not material_mass:
            m.add_bar_mass(b, kind="uniform", m1=rho * A)      # translation only: Euler-Bernoulli
    for x in nd:
        m.add_support(x, dx=True, dy=True, rx=True, rz=True)   # bending in the XZ plane
    m.add_support(nd[0], dz=True)
    m.add_support(nd[-1], dz=True)
    middle = nd[n // 2]
    if attachment == "spring":
        m.add_support(middle, kz=k_spring)
    elif attachment == "support":
        m.add_support(middle, dz=True)
    elif attachment == "mass":
        m.add_point_mass(middle, mass)
    return m.modal(n_modes=3, method="dense").frequencies

def f_of(bl):                                # frequency associated with beta l, l = L / 2
    return (2 * bl / L) ** 2 / (2 * math.pi) * math.sqrt(E * I / (rho * A))

def g(bl, attachment):                       # symmetric modes (half-beam 0 <= x <= l)
    s, c, sh, ch = math.sin(bl), math.cos(bl), math.sinh(bl), math.cosh(bl)
    if attachment == "support":
        return s * ch - c * sh               # tan(beta l) = tanh(beta l)
    kr = {"spring": k_spring, "mass": -mass * (2 * math.pi * f_of(bl)) ** 2}.get(attachment, 0.0)
    return 4 * E * I * (2 * bl / L) ** 3 * c * ch + kr * (s * ch - c * sh)

def exact(attachment):
    f, x, step = [], 0.1, 1e-3
    while x < 3 * math.pi:                   # scan in beta l, then bisection on each root
        if g(x, attachment) * g(x + step, attachment) < 0:
            u, v = x, x + step
            for _ in range(60):
                if g(u, attachment) * g((u + v) / 2, attachment) <= 0:
                    v = (u + v) / 2
                else:
                    u = (u + v) / 2
            f.append(f_of(u))
        x += step
    f += [f_of(j * math.pi / 2) for j in (2, 4)]      # antisymmetric: beta L = 2 pi, 4 pi
    return sorted(f)[:3]

def row(name, fk, fe, ref):
    diff = max(abs(u / v - 1) for u, v in zip(fk, fe))
    print(f"{name:8s} keta     " + "".join(f"{v:11.4f}" for v in fk))
    print(f"{'':8s} {ref:9s}" + "".join(f"{v:11.4f}" for v in fe) + f"  diff. {100 * diff:.4f} %")
    return diff

print(f"{'attachment':18s}" + "".join(f"{'f%d [Hz]' % j:>11s}" for j in (1, 2, 3)))
for attachment in ("none", "spring", "support", "mass"):
    assert row(attachment, keta_f(attachment), exact(attachment), "exact") < 1e-5

# Material mass: the rotary inertia of the sections is included (Rayleigh beam)
fr = [f_of(j * math.pi / 2) / math.sqrt(1 + (j * math.pi / L) ** 2 * I / A) for j in (1, 2, 3)]
assert row("rho", keta_f("none", material_mass=True), fr, "Rayleigh") < 5e-5
