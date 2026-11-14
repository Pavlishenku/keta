# COMPLETE EXAMPLE 1 - Closed reinforced concrete box frame on elastic soil (hydraulic structure
# or underpass), plane analysis in XZ, slice of 1 m.
#   span between axes 6.0 m, height between axes 4.0 m, thicknesses 0.50 m;
#   subgrade reaction modulus of the soil 30 MN/m3; fill of 1.0 m on the top slab.
import keta

B, H, e, k_soil = 6.0, 4.0, 0.50, 30e6
gb, gs, K0 = 25e3, 20e3, 0.5                     # unit weights of concrete / fill, earth pressure coefficient
nb = 12                                          # bars per slab, 8 per side wall

m = keta.Model()
concrete = m.add_material(e=32e9, nu=0.2, rho=gb / 9.80665)
sec = m.add_section_rect(b=1.0, h=e)
bot = [m.add_node(B * i / nb, 0.0, 0.0) for i in range(nb + 1)]           # base slab
top = [m.add_node(B * i / nb, 0.0, H) for i in range(nb + 1)]             # top slab
left = [bot[0]] + [m.add_node(0.0, 0.0, H * i / 8) for i in range(1, 8)] + [top[0]]
right = [bot[-1]] + [m.add_node(B, 0.0, H * i / 8) for i in range(1, 8)] + [top[-1]]
base = [m.add_beam(bot[i], bot[i + 1], sec, concrete) for i in range(nb)]
roof = [m.add_beam(top[i], top[i + 1], sec, concrete) for i in range(nb)]
wl = [m.add_beam(left[i], left[i + 1], sec, concrete) for i in range(8)]
wr = [m.add_beam(right[i], right[i + 1], sec, concrete) for i in range(8)]
for b in base:
    m.set_soil(b, kfz=k_soil * 1.0)                                        # soil under the base slab
for x in set(bot + top + left + right):
    m.add_support(x, dy=True, rx=True, rz=True)                            # plane analysis XZ
m.add_support(bot[nb // 2], dx=True)                                       # horizontal restraint

# Case 1: self-weight; case 2: fill on the top slab; case 3: earth pressure (at rest);
# case 4: live load of 10 kN/m2 on the top slab and the corresponding earth pressure
m.add_self_weight(case=1)
for b in roof:
    m.add_udl(b, wz=-gs * 1.0, case=2)
for i in range(8):
    z1, z2 = H * i / 8, H * (i + 1) / 8
    p1, p2 = K0 * gs * (1.0 + H - z1), K0 * gs * (1.0 + H - z2)            # pressure at depth
    m.add_linear_load(wl[i], w1=(p1, 0, 0), w2=(p2, 0, 0), case=3)
    m.add_linear_load(wr[i], w1=(-p1, 0, 0), w2=(-p2, 0, 0), case=3)
    m.add_udl(wl[i], wx=K0 * 10e3, case=4)
    m.add_udl(wr[i], wx=-K0 * 10e3, case=4)
for b in roof:
    m.add_udl(b, wz=-10e3, case=4)
names = ["self-weight", "fill", "earth pressure", "live load"]
for c, name in enumerate(names, start=1):
    m.name_case(c, name)
uls = m.add_combination([(1.35, 1), (1.35, 2), (1.35, 3), (1.5, 4)], name="ULS")
sls = m.add_combination([(1.0, 1), (1.0, 2), (1.0, 3), (1.0, 4)], name="characteristic SLS")

r = m.analyze(stations=5)
print(m)
print("problems:", m.validate())
print(f"{'':22} {'M top corner':>13} {'M mid-top':>11} {'M bot corner':>13} {'M mid-base':>12} {'N wall':>11}")
kN = lambda v: round(v / 1e3, 2) + 0.0                # kN or kN.m, without negative zero
for name, res in [(names[c - 1], r.case(c)) for c in r.case_ids] + [("ULS", r.combo(uls)), ("SLS", r.combo(sls))]:
    Ms = res.bar_forces(roof[0])[0][5]
    Mt = res.bar_forces(roof[nb // 2])[0][5]
    Mi = res.bar_forces(base[0])[0][5]
    Mr = res.bar_forces(base[nb // 2])[0][5]
    N = res.bar_forces(wl[0])[0][1]
    print(f"{name:22} {kN(Ms):13.2f} {kN(Mt):11.2f} {kN(Mi):13.2f} {kN(Mr):12.2f} {kN(N):11.2f}")

# Soil pressure under the base slab at SLS: p = k w (computed from the displacements)
res = r.combo(sls)
p = [-k_soil * res.displacement(x)[2] / 1e3 for x in bot]
print("soil pressure under the base slab, SLS [kPa]:", " ".join(f"{v:.1f}" for v in p[::2]))
print(f"mean settlement {sum(res.displacement(x)[2] for x in bot) / len(bot) * 1e3:.3f} mm; "
      f"resultant of the pressures {sum((p[i] + p[i + 1]) / 2 * B / nb for i in range(nb)):.2f} kN "
      f"(weight + loads: {(gb * e * (2 * B + 2 * H) + gs * B + 10e3 * B) / 1e3:.2f} kN)")
