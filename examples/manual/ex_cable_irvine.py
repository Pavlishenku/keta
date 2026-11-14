# Heavy cable of 100 m (40 catenaries, 10 kg/m, EA = 1e9 N), sag close to L/100: modes
# about the equilibrium (tangent of the large-displacement state) against Irvine's theory.
import math
import keta

L, N, MU, EA = 100.0, 40, 10.0, 1e9
d0 = L / 100
H0 = MU * 9.80665 * L * L / (8 * d0)
l0 = L / N * (1 + 8 / 3 * (d0 / L) ** 2 - H0 / EA)        # unstressed length of one element
m = keta.Model()
nd = [m.add_node(L * i / N, 0, 0) for i in range(N + 1)]
for i in range(N):
    m.add_cable(nd[i], nd[i + 1], ea=EA, w=MU * 9.80665, l0=l0)
m.fix(nd[0])
m.fix(nd[N])
for x in nd[1:N]:
    m.add_support(x, rx=True, ry=True, rz=True)
m.name_case(1, "self-weight")                            # empty case: the weight of the cables alone

state = m.nonlinear(case=1, substeps=20, max_iterations=80)
sag = -state.displacement(nd[N // 2])[2]
modes = state.tangent().modal(n_modes=6)
H = MU * 9.80665 * L * L / (8 * sag)                     # horizontal tension (parabola)
mu = MU * N * l0 / L                                      # mass per metre of span
scale = L / math.sqrt(H / mu)
print(f"sag {sag:.4f} m, H = {H / 1e3:.2f} kN")
for i in range(modes.n_modes):
    dy = max(abs(d[1]) for _, d in modes.mode(i))
    dz = max(abs(d[2]) for _, d in modes.mode(i))
    print(f"mode {i + 1}: f = {modes.frequencies[i]:.4f} Hz, w L / sqrt(H/m) = {modes.omegas[i] * scale:.3f}, "
          f"{'out-of-plane' if dy > dz else 'in-plane'}")
lam2 = (mu * 9.80665 * L / H) ** 2 * L / (H * L * (1 + 8 * (sag / L) ** 2) / EA)
lo, hi = math.pi / 2 + 1e-12, 3 * math.pi / 2 - 1e-12          # tan x = x - (4 / lambda2) x^3
for _ in range(200):
    x = 0.5 * (lo + hi)
    lo, hi = (x, hi) if math.tan(x) - x + 4 / lam2 * x ** 3 < 0 else (lo, x)
print(f"lambda2 = {lam2:.1f}; Irvine: out-of-plane n pi = 3.142, 6.283, 9.425; antisymmetric 2 pi = 6.283; "
      f"symmetric {lo + hi:.3f}")
