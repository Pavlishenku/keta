# Spectre de l'EN 1998-1 et modes jusqu'à une fréquence de coupure.
import math
import keta

ag, S, TB, TC, TD = 2.0, 1.2, 0.15, 0.5, 2.0
spectre = keta.ec8_spectrum(ag, S, TB, TC, TD)                 # élastique, 5 %
calcul = keta.ec8_spectrum(ag, S, TB, TC, TD, q=1.5)           # de calcul
eta = math.sqrt(10 / (5 + 100 * 0.05))
print(f"plateau élastique = {dict(spectre)[0.3]:.3f} m/s2 (ag S eta 2,5 = {ag * S * eta * 2.5:.3f})")
print(f"plateau de calcul  = {dict(calcul)[0.3]:.3f} m/s2 (ag S 2,5/q = {ag * S * 2.5 / 1.5:.3f})")

# pile de 12 m encastrée en pied, masse de tablier en tête
m = keta.Model()
n = [m.add_node(0, 0, 12.0 * i / 6) for i in range(7)]
s = m.add_section_rect(b=2.0, h=2.0)
for i in range(6):
    m.add_beam(n[i], n[i + 1], s, m.concrete())
m.fix(n[0])
m.add_point_mass(n[-1], 400e3)
mo = m.modal(f_max=15.0)
print(f"{len(mo.frequencies)} modes sous 15 Hz :", [round(f, 3) for f in mo.frequencies])
sp = m.response_spectrum(spectre, dir=0, zeta=0.05, rule="cqc", n_modes=len(mo.frequencies))
print(f"effort tranchant à la base (X) = {sp.base_shear() / 1e3:.1f} kN")
