# L-shaped frame fixed at A, on a sliding support at C: outer faces at 0 °C,
# inner faces at +10 °C.
import keta

a, b, h = 4.0, 0.3, 0.4                   # length of the post and of the girder, section b x h
E, alpha, t_out, t_in = 30e9, 1e-5, 0.0, 10.0
I = b * h**3 / 12

m = keta.Model()
concrete = m.add_material(e=E, nu=0.2, alpha_temp=alpha)
sec = m.add_section_props(sx=1e3, iy=I, iz=I, ix=I)      # axially rigid, as in the hand calculation
A, B, C = m.add_node(0, 0, 0), m.add_node(0, 0, a), m.add_node(a, 0, a)
post, girder = m.add_beam(A, B, sec, concrete), m.add_beam(B, C, sec, concrete)
m.fix(A)
m.add_support(C, dz=True)
for n in (B, C):
    m.add_support(n, dy=True, rx=True, rz=True)          # plane frame in XZ
# local z: -X for the post, +Z for the girder; the inner face is on the z < 0 side of both bars
for bar in (post, girder):
    m.add_thermal(bar, dt=(t_out + t_in) / 2, gz=-(t_in - t_out) / h)
case = m.analyze().case(1)

X1 = -dict(case.reactions())[C][2]                    # reaction at C, positive downwards
MB = case.bar_forces(girder)[0][5]                    # My > 0: top (cold) fibre in tension
sigma = MB * (h / 2) / I
X1_ref = 15 * alpha * E * I / (4 * a**2) * (3 * a / h + 1)
for name, v, ref in (("reaction X1 at C [N]", X1, X1_ref),
                     ("moment at B [N.m]", MB, a * X1_ref),
                     ("stress in cold fibre at B [kPa]", sigma / 1e3, a * X1_ref * h / 2e3 / I)):
    diff = 100 * (v - ref) / ref
    print(f"{name:32s} keta {v:10.4f}  closed form {ref:10.4f}  diff. {round(diff, 4) + 0:.4f} %")
    assert abs(diff) < 1e-3
