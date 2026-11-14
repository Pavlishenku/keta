# Staged construction: two precast spans of 20 m placed as simply supported, then made continuous
# by a closure over the central support; creep redistribution of the support moment.
import math
import keta

L, q, E = 20.0, -50e3, 35e9
m = keta.Model()
concrete = m.add_material(e=E, nu=0.2, rho=0.0)
m.set_concrete_law(concrete, "COURBON", e_direct=E)   # non-ageing creep, phi(inf) = 2
sec = m.add_section_props(sx=0.8, iy=0.25, iz=0.1, ix=0.1)   # rigid in shear
n = [m.add_node(L * i / 4, 0, 0) for i in range(9)]
b = [m.add_beam(n[i], n[i + 1], sec, concrete) for i in range(8)]
m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(n[4], dy=True, dz=True)
m.add_support(n[8], dy=True, dz=True)
for bi in b:
    m.add_udl(bi, wz=q, case=1)                    # self-weight + superimposed dead loads

phases = [
    dict(name="placing", time=20, place_hinges=[dict(bar=b[3], end=1, dofs="ry")]),
    dict(name="loading", time=28, cases=[1]),
    dict(name="closure", time=60, close_hinges=[dict(bar=b[3], end=1, dofs="ry")]),
    dict(name="1 year", time=393),
    dict(name="30 years", time=28 + 30 * 365),
]
r = m.phased(phases, creep_along_bar=True, creep_substeps_per_decade=32)
print("stages:", r.names())
for i in range(len(r)):
    p = r.phase(i)
    Ma = p.bar_forces(b[3])[-1][5] + 0.0         # moment at the end of the bar over the support
    print(f"{p.name:11} t = {p.time:6.0f} days: M(support) = {abs(Ma) / 1e3:8.2f} kN.m, "
          f"midspan deflection = {p.displacement(n[2])[2] * 1e3:8.3f} mm")
r0 = m.phased(phases)                            # default options: creep at mid-bar
print(f"default options, 30 years: M(support) = {r0.phase(4).bar_forces(b[3])[-1][5] / 1e3:.2f} kN.m")

# Check: non-ageing law, exact limit M = Mc (phi(inf) - phi(t1 - t0)) / (1 + phi(inf))
phi = lambda d: 2 * (1 - math.exp(-0.006 * d))
Mc = -q * L * L / 8                              # moment of the beam cast continuous
M30 = Mc * (phi(30 * 365) - phi(60 - 28)) / (1 + phi(30 * 365))
print(f"beam cast continuous: Mc = {Mc / 1e3:.2f} kN.m; limit by creep: {M30 / 1e3:.2f} kN.m")
