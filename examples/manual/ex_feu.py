# Fire analysis: fire curves, properties of heated materials, temperature
# of a steel profile, then a 6 m RC beam exposed to ISO 834 on its soffit and sides.
import math
import keta

for curve in ("iso834", "external", "hydrocarbon"):
    print(f"{curve:12}: " + ", ".join(f"{t} min {keta.fire_gas_temperature(curve, t):6.1f} °C"
                                      for t in (15, 30, 60, 120)))
print(f"steel at 500 °C: ky = {keta.steel_ky(500):.2f}, kE = {keta.steel_ke(500):.2f}; "
      f"critical temperature (mu0 = 0.5) = {keta.steel_critical_temperature(0.5):.1f} °C")
print(f"siliceous concrete at 500 °C: kc = {keta.concrete_kc(500):.2f}, eps_c1 = {keta.concrete_eps_c1(500):.4f}, "
      f"thermal strain = {keta.concrete_thermal_strain(500):.5f}")
print(f"reinforcement (class N, cold-worked) at 500 °C: ks = {keta.rebar_ks(500):.3f}, kE = {keta.rebar_kE(500):.3f}")
T = keta.steel_temperature_unprotected("iso834", am_v=200.0, ksh=1.0, t_end_minutes=30.0, dt=5.0)
print(f"unprotected profile Am/V = 200 m-1: {T[-1]:.1f} °C after 30 min ({len(T)} steps of 5 s)")

def make_beam(q):
    """Beam 0.20 x 0.30 m over 6 m, 3 HA16 at the top and at the bottom; beta = pi/2 makes
    the thickness of the fire section (local y axis) vertical, face 'sup' at the top."""
    m = keta.Model()
    concrete = m.add_material(e=30e9, nu=0.2, rho=0.0)
    s = m.add_section_rect(0.2, 0.3)
    n = [m.add_node(1.0 * i, 0, 0) for i in range(7)]
    b = [m.add_beam(n[i], n[i + 1], s, concrete, beta=math.pi / 2) for i in range(6)]
    m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
    m.add_support(n[-1], dy=True, dz=True)
    for bi in b:
        m.add_udl(bi, wz=-q, case=1)
    layer = lambda face: dict(diameter=0.016, count=3, placement=face, cover=0.04)
    fire_bars = [dict(bar=bi, fck=30e6, rebars=[layer("inf"), layer("sup")],
                      section=dict(kind="rectangle", thickness=0.3, width=0.2,
                                   exposures=dict(inf=0, pz=0, mz=0), mesh=0.02)) for bi in b]
    return m, n, fire_bars

m, n, fire_bars = make_beam(10e3)
cold = m.analyze().case(1).displacement(n[3])[2]
for duration in (30, 45, 60):
    try:
        r = m.fire(fire_bars, [dict(kind="iso834")], duration=duration, preload_case=1, local=True)
    except ValueError:
        print(f"ISO 834, {duration} min: no equilibrium, failure of the beam before {duration} min")
        continue
    s = r.sections()[2]
    print(f"ISO 834, {duration} min: deflection {r.displacement(n[3])[2] * 1e3:8.2f} mm (cold {cold * 1e3:.2f}), "
          f"{r.iterations} iterations, sum of the reactions {sum(x[1][2] for x in r.reactions()) / 1e3:.3f} kN")
    print(f"   section: theta max {s['theta_max']:.0f} °C, mean {s['theta_mean']:.0f} °C, "
          f"reinforcement {s['theta_rebar_max']:.0f} °C")
    print(f"   hot stiffnesses: EA = {s['EA'] / 1e6:.1f} MN, EIy = {s['EIy'] / 1e6:.3f} MN.m2, "
          f"EIz = {s['EIz'] / 1e6:.3f} MN.m2, centre (yc, zc) = ({round(s['yc'] * 1e3, 1) + 0.0:.1f}, {round(s['zc'] * 1e3, 1) + 0.0:.1f}) mm")

# Failure: after 1 min of fire the section is at 20 °C; Mu = 70.8 kN.m, mean moment of the central
# bars 13 q / 3 -> q_u = 3 Mu / 13
q_u = 3 * 70.8e3 / 13
for f in (0.97, 1.03):
    m, n, fire_bars = make_beam(f * q_u)
    try:
        m.fire(fire_bars, [dict(kind="iso834")], duration=1, preload_case=1, local=True)
        print(f"q = {f:.2f} q_u: equilibrium found")
    except ValueError:
        print(f"q = {f:.2f} q_u: no equilibrium (failure)")
