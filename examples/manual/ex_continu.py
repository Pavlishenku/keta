# Continuum elements: beam modelled as a plane-stress wall (Q8 quadrilaterals),
# clamped, with a distributed shear force of 20 kN at its free end.
import keta

L, h, t, P, E = 4.0, 0.4, 0.1, 20e3, 30e9
nx, ny = 20, 2
m = keta.Model()
mat = m.add_material(e=E, nu=0.0)
# nodes at half-element spacing (the element centres are not used)
N = {}
for j in range(2 * ny + 1):
    for i in range(2 * nx + 1):
        if i % 2 == 0 or j % 2 == 0:
            N[i, j] = m.add_node(L * i / (2 * nx), h * j / (2 * ny), 0.0)
elements = []
for j in range(0, 2 * ny, 2):
    for i in range(0, 2 * nx, 2):
        corners = [N[i, j], N[i + 2, j], N[i + 2, j + 2], N[i, j + 2]]            # anticlockwise
        midsides = [N[i + 1, j], N[i + 2, j + 1], N[i + 1, j + 2], N[i, j + 1]]   # sides 01, 12, 23, 30
        elements.append(m.add_continuum("plane_stress", corners + midsides, mat, t=t))
for j in range(2 * ny + 1):
    m.add_support(N[0, j], dx=True, dy=True)                   # clamping of the end face x = 0
for e in elements[nx - 1::nx]:                                 # last element of each row
    m.add_continuum_traction(e, face=1, traction=(0.0, -P / (h * t), 0.0))   # side x = L

case_ = m.analyze().case(1)
v = case_.displacement(N[2 * nx, ny])[1]
I = t * h**3 / 12
print(f"tip deflection (mid-depth fibre): {v * 1e3:.4f} mm")
print(f"Timoshenko beam (shear 5/6): {-(P * L**3 / (3 * E * I) + P * L / (5 / 6 * E / 2 * t * h)) * 1e3:.4f} mm")
sx, sy, txy, sz = case_.continuum_stress(elements[nx // 2], (0.0, -1.0, 0.0))  # midspan, bottom edge
print(f"sigma_x at x = {L / 2 + L / (2 * nx):.1f} m, bottom fibre: {sx / 1e6:.3f} MPa "
      f"(simple bending: {-P * (L / 2 - L / (2 * nx)) * h / 2 / I / 1e6:.3f} MPa)")
