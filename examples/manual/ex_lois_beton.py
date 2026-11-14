# Time-dependent concrete laws: shortening of a 4 m column under 2 MN applied
# at 28 days, for several creep and shrinkage laws.
import math
import keta

def column(rule=None, functions=None, loaded=True, **params):
    m = keta.Model()
    c = m.add_material(e=35e9, nu=0.2, rho=0.0)
    if rule:
        m.set_concrete_law(c, rule, **params)
    if functions:
        m.set_concrete_law_functions(c, *functions)
    s = m.add_section_rect(0.4, 0.4)
    a, b = m.add_node(0, 0, 0), m.add_node(0, 0, 4)
    m.add_beam(a, b, section=s, material=c)
    m.fix(a)
    m.add_support(b, dx=True, dy=True, rx=True, ry=True, rz=True)
    m.add_point_load(b, fz=-2e6 if loaded else 0.0, case=1)
    dates = [28, 128, 365, 3650, 10000]
    r = m.phased([dict(name=f"{t} days", time=t, cases=[1] if t == 28 else []) for t in dates],
                 creep_substeps_per_decade=16)
    return [r.phase(i).displacement(b)[2] for i in range(len(r))]

u0 = -2e6 * 4 / (35e9 * 0.16)
print(f"elastic shortening: {u0 * 1e3:.4f} mm")
print("dates [days]          :       28      128      365     3650    10000")
cases = [("COURBON", dict(e_direct=35e9)),
         ("EC2", dict(fc=35e6, rh=70.0, h_mm=400.0, cement="N")),
         ("EC2, creep only", None),
         ("BPEL", dict(fc=35e6, rh=70.0, rm_cm=10.0)),
         ("FIB2010", dict(fc=35e6, rh=70.0, h_mm=400.0, cement="N"))]
for name, kw in cases:
    if kw is None:    # shrinkage alone (unloaded column) subtracted from the EC2 law
        u = [a - b for a, b in zip(column("EC2", **cases[1][1]),
                                   column("EC2", loaded=False, **cases[1][1]))]
    else:
        u = column(name, **kw)
    print(f"{name:17} [mm]: " + " ".join(f"{x * 1e3:8.3f}" for x in u))

# User-defined law: phi(t, t0) = 1.5 (1 - exp(-(t - t0) / 100)), constant modulus
phi = lambda t, t0: 1.5 * (1 - math.exp(-(t - t0) / 100)) if t > t0 else 0.0
u = column(functions=(phi, lambda t0: 35e9))
print(f"{'user-defined':17} [mm]: " + " ".join(f"{x * 1e3:8.3f}" for x in u))
print(f"COURBON check at 10000 days: u0 (1 + phi) = {u0 * (1 + 2 * (1 - math.exp(-0.006 * 9972))) * 1e3:.3f} mm")
