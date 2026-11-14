# Charges réparties sur une barre inclinée : repère global, local et projeté.
import math
import keta

alpha = math.radians(30)
L = 10.0
for repere in ("global", "local", "projected"):
    m = keta.Model()
    a = m.add_node(0, 0, 0)
    b = m.add_node(L * math.cos(alpha), 0, L * math.sin(alpha))
    barre = m.add_beam(a, b, m.add_section_rect(0.3, 0.6), m.concrete())
    m.add_support(a, dx=True, dy=True, dz=True, rx=True)
    m.add_support(b, dy=True, dz=True)
    m.add_udl(barre, wz=-10e3, frame=repere)          # 10 kN/m selon -Z (ou -z local)
    Fx, Fy, Fz = m.analyze().case(1).sum_reactions()
    print(f"{repere:9} : résultante des réactions Fx = {round(Fx / 1e3, 3) + 0.0:7.3f} kN, Fz = {Fz / 1e3:7.3f} kN")
print(f"attendu : global 10 L = {10 * L:.3f} kN ; projeté 10 L cos(alpha) = {10 * L * math.cos(alpha):.3f} kN ; "
      f"local : 10 L normal à la barre")
