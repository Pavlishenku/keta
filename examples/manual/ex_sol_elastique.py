# Beam on elastic soil (Winkler): ground beam of 40 m under a centred load of 500 kN,
# compared with Hetényi's infinite beam: the solution is exact whatever the mesh.
import keta

L, P, k = 40.0, 500e3, 40e6               # length [m], load [N], modulus x width [N/m per m]
EI = 30e9 * 1.0 * 0.8**3 / 12
lam = (k / (4 * EI)) ** 0.25
print(f"characteristic length 1/lambda = {1 / lam:.3f} m")
print(f"infinite beam: w0 = {-P * lam / (2 * k) * 1e3:.4f} mm, M0 = {P / (4 * lam) / 1e3:.3f} kN.m")

for n in (2, 8, 40):
    m = keta.Model()
    nd = [m.add_node(L * i / n, 0, 0) for i in range(n + 1)]
    sec = m.add_section_props(sx=0.8, iy=1.0 * 0.8**3 / 12, iz=0.8 * 1.0**3 / 12, ix=0.1)  # rigid in shear
    mat = m.add_material(e=30e9, nu=0.2)
    bars = [m.add_beam(nd[i], nd[i + 1], sec, mat) for i in range(n)]
    for b in bars:
        m.set_soil(b, kfz=k)              # vertical elastic foundation (local z = Z)
    mid = nd[n // 2]
    m.add_support(mid, dx=True, dy=True, rx=True, rz=True)      # out-of-plane movements restrained
    m.add_point_load(mid, fz=-P)
    case_ = m.analyze().case(1)
    w0 = case_.displacement(mid)[2]
    M0 = case_.bar_forces(bars[n // 2])[0][5]
    assert abs(case_.sum_reactions()[2]) < 1e-3     # the whole load goes into the soil
    print(f"{n:4d} bars: w0 = {w0 * 1e3:.4f} mm, M0 = {-M0 / 1e3:.3f} kN.m")
