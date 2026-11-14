# Bar with variable section: cantilever of 8 m with linearly varying depth (1.0 -> 0.5 m),
# then inverse parabolic second moment of area; comparison with the integral of M m / EI(x).
import keta

E, L, P, b = 30e9, 8.0, 100e3, 0.3

def cantilever(law, orig, ext, mid=None, n=1):
    m = keta.Model()
    c = m.add_material(e=E, nu=0.2, rho=0.0)
    s = m.add_section_props(sx=0.3, iy=0.025, iz=0.002, ix=0.005)    # "mean" section of the bar
    nd = [m.add_node(L * i / n, 0, 0) for i in range(n + 1)]
    bars = [m.add_beam(nd[i], nd[i + 1], s, c) for i in range(n)]
    m.fix(nd[0])
    m.set_tapered(bars[0], law, orig, ext, mid)
    m.add_point_load(nd[-1], fz=-P)
    return m.analyze().case(1).displacement(nd[-1])[2]

def rect(h):                                    # properties of a b x h rectangle (no shear)
    return dict(A=b * h, Iy=b * h**3 / 12, Iz=h * b**3 / 12, J=0.1 * b**3 * h)

def simpson(f, a, c, n=2000):
    step = (c - a) / n
    return step / 3 * sum((1 if i in (0, n) else 4 if i % 2 else 2) * f(a + i * step) for i in range(n + 1))

w = cantilever("lin_z", rect(1.0), rect(0.5))
h = lambda x: 1.0 + (0.5 - 1.0) * x / L
wref = simpson(lambda x: P * (L - x) ** 2 / (E * b * h(x) ** 3 / 12), 0, L)
print(f"linear depth: uz = {w * 1e3:.5f} mm; integral of M m / EI = {-wref * 1e3:.5f} mm")

I0, Im, I1 = 0.025, 0.012, 0.004
w = cantilever("para", dict(A=0.3, Iy=I0, Iz=0.002, J=0.005), dict(A=0.2, Iy=I1, Iz=0.002, J=0.005),
               dict(A=0.25, Iy=Im, Iz=0.002, J=0.005))
inv = lambda t: (1 / I0) * 2 * (t - 0.5) * (t - 1) - (1 / Im) * 4 * t * (t - 1) + (1 / I1) * 2 * t * (t - 0.5)
wref = simpson(lambda x: P * (L - x) ** 2 * inv(x / L) / E, 0, L)
print(f"inverse parabolic: uz = {w * 1e3:.5f} mm; integral = {-wref * 1e3:.5f} mm")

# the same cantilever split into 8 prismatic bars of mean depth
m = keta.Model()
c = m.add_material(e=E, nu=0.2, rho=0.0)
nd = [m.add_node(i * 1.0, 0, 0) for i in range(9)]
for i in range(8):
    hm = h(i + 0.5)
    m.add_beam(nd[i], nd[i + 1], m.add_section_props(sx=b * hm, iy=b * hm**3 / 12, iz=hm * b**3 / 12, ix=0.01), c)
m.fix(nd[0])
m.add_point_load(nd[-1], fz=-P)
print(f"8 prismatic bars: uz = {m.analyze().case(1).displacement(nd[-1])[2] * 1e3:.5f} mm")
