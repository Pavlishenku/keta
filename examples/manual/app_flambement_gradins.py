# Buckling of a telescopic column fixed at its foot, free at its top, made of 1 to 8 segments of
# decreasing section, compressed by a vertical force at the top.
import math
import keta

E = 206e9
SEGMENTS = [  # (length [m], area [m2], inertia 1 [m4], inertia 2 [m4]), from the foot to the top
    (15.3, 0.0906, 0.04326, 0.03466), (15.5, 0.0645, 0.02609, 0.02281),
    (15.4, 0.0529, 0.01906, 0.01703), (15.6, 0.0478, 0.01441, 0.01335),
    (15.4, 0.0415, 0.01039, 0.01043), (15.3, 0.0346, 0.00749, 0.00756),
    (15.1, 0.0296, 0.00541, 0.00560), (14.8, 0.0302, 0.00463, 0.00489)]

def keta_pcr(n, subdivision=4):
    m = keta.Model()
    steel = m.add_material(e=E, nu=0.3)
    foot = m.add_node(0, 0, 0)
    m.fix(foot)
    z, prev = 0.0, foot
    for L, A, I1, I2 in SEGMENTS[:n]:
        sec = m.add_section_props(sx=A, iy=I1, iz=I2, ix=I1 + I2)
        for _ in range(subdivision):                # 4 bars per segment
            z += L / subdivision
            nn = m.add_node(0, 0, z)
            m.add_beam(prev, nn, sec, steel)
            prev = nn
    m.add_point_load(prev, fz=-1e3)                 # 1 kN: the critical factor is Pcr in kN
    return m.buckling(n_modes=1).factors[0]

def residual(P, n, plane):
    # u = top offset - deflection: u'' + (P / EI) u = 0 in each segment,
    # u = 1 and u' = 0 at the foot (fixed), u = 0 at the top (point of application of P)
    u, du = 1.0, 0.0
    for t in SEGMENTS[:n]:
        L, k = t[0], math.sqrt(P / (E * t[2 + plane]))
        c, s = math.cos(k * L), math.sin(k * L)
        u, du = u * c + du / k * s, -u * k * s + du * c
    return u

def pcr_exact(n):
    values = []
    for plane in (0, 1):                            # each of the two principal planes
        P = 1e3
        while residual(P * 1.01, n, plane) > 0:     # first root: first change of sign
            P *= 1.01
        a, b = P, P * 1.01                          # then bisection
        for _ in range(60):
            c = (a + b) / 2
            if residual(c, n, plane) > 0:
                a = c
            else:
                b = c
        values.append(a)
    return min(values) / 1e3

L1, I2 = SEGMENTS[0][0], SEGMENTS[0][3]                # smaller inertia of segment 1
print(f"1 segment, Euler pi2 E I / 4 L2 = {math.pi**2 * E * I2 / (4 * L1**2) / 1e3:.1f} kN")
print(" n    height     Pcr keta       Pcr exact     diff.")
for n in range(1, 9):
    pk, pe = keta_pcr(n), pcr_exact(n)
    diff = 100 * (pk - pe) / pe
    H = sum(t[0] for t in SEGMENTS[:n])
    print(f"{n:2d}   {H:6.1f} m  {pk:10.1f} kN  {pe:10.1f} kN  {round(diff, 4) + 0:7.4f} %")
    assert abs(diff) < 0.01
