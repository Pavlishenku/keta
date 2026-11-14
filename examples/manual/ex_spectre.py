# Response spectrum analysis: 3-storey frame with rigid floors, horizontal design
# spectrum of the EN 1998-1 type (built in Python), SRSS and CQC combinations.
import keta

def design_spectrum(ag, S=1.2, TB=0.15, TC=0.5, TD=2.0, q=1.5, beta=0.2):
    """Design spectrum (m/s2) tabulated between 0 and 4 s."""
    pts = []
    for i in range(401):
        T = 0.01 * i
        if T <= TB:
            sa = ag * S * (2 / 3 + T / TB * (2.5 / q - 2 / 3))
        elif T <= TC:
            sa = ag * S * 2.5 / q
        elif T <= TD:
            sa = max(ag * S * 2.5 / q * TC / T, beta * ag)
        else:
            sa = max(ag * S * 2.5 / q * TC * TD / T**2, beta * ag)
        pts.append((T, sa))
    return pts

m = keta.Model()
concrete = m.add_material(e=30e9, nu=0.2, rho=0.0)             # mass carried by the floors
col = m.add_section_rect(0.4, 0.4)
pts = [[m.add_node(x, 0, 3.0 * k) for x in (0.0, 6.0)] for k in range(4)]
for k in range(3):
    for j in range(2):
        m.add_beam(pts[k][j], pts[k + 1][j], col, concrete)
for j in range(2):
    m.fix(pts[0][j])
for k in range(1, 4):
    m.add_rigid_link(pts[k][0], [pts[k][1]], dy=False, dz=False, rx=False, ry=False, rz=False)  # floor
    for j in range(2):
        m.add_support(pts[k][j], dy=True, rx=True, ry=True, rz=True)   # plane frame XZ, infinitely stiff beams
        m.add_point_mass(pts[k][j], 30e3)

mo = m.modal(n_modes=3)
print("periods [s]:", [round(T, 4) for T in mo.periods])
sp = design_spectrum(ag=2.0)
for rule in ("srss", "cqc"):
    r = m.response_spectrum(spectrum=sp, dir=0, rule=rule, zeta=0.05, n_modes=3)
    print(f"{rule.upper():4}: base shear {r.base_shear() / 1e3:8.2f} kN, "
          f"top displacement {r.combined_disp(pts[3][0])[0] * 1e3:7.3f} mm")
print("modal base shears [kN]:", [round(v / 1e3, 2) for v in r.modal_base_shear()])
print("effective masses [t]:", [round(x / 1e3, 3) for x in r.effective_mass()], "/",
      round(r.total_mass() / 1e3, 3))
