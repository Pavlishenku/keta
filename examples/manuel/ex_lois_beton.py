# Lois différées du béton : raccourcissement d'un poteau de 4 m sous 2 MN appliqués
# à 28 jours, pour plusieurs lois de fluage et de retrait.
import math
import keta

def poteau(regle=None, fonctions=None, charge=True, **parametres):
    m = keta.Model()
    c = m.add_material(e=35e9, nu=0.2, rho=0.0)
    if regle:
        m.set_concrete_law(c, regle, **parametres)
    if fonctions:
        m.set_concrete_law_functions(c, *fonctions)
    s = m.add_section_rect(0.4, 0.4)
    a, b = m.add_node(0, 0, 0), m.add_node(0, 0, 4)
    m.add_beam(a, b, section=s, material=c)
    m.fix(a)
    m.add_support(b, dx=True, dy=True, rx=True, ry=True, rz=True)
    m.add_point_load(b, fz=-2e6 if charge else 0.0, case=1)
    dates = [28, 128, 365, 3650, 10000]
    r = m.phased([dict(name=f"{t} j", time=t, cases=[1] if t == 28 else []) for t in dates],
                 creep_substeps_per_decade=16)
    return [r.phase(i).displacement(b)[2] for i in range(len(r))]

u0 = -2e6 * 4 / (35e9 * 0.16)
print(f"raccourcissement élastique : {u0 * 1e3:.4f} mm")
print("dates [j]                 :      28      128      365     3650    10000")
cas = [("COURBON", dict(e_direct=35e9)),
       ("EC2", dict(fc=35e6, rh=70.0, h_mm=400.0, cement="N")),
       ("EC2, fluage seul", None),
       ("BPEL", dict(fc=35e6, rh=70.0, rm_cm=10.0)),
       ("FIB2010", dict(fc=35e6, rh=70.0, h_mm=400.0, cement="N"))]
for nom, par in cas:
    if par is None:   # retrait seul (poteau non chargé) retranché de la loi EC2
        u = [a - b for a, b in zip(poteau("EC2", **cas[1][1]),
                                   poteau("EC2", charge=False, **cas[1][1]))]
    else:
        u = poteau(nom, **par)
    print(f"{nom:17} [mm] : " + " ".join(f"{x * 1e3:8.3f}" for x in u))

# Loi définie par l'utilisateur : phi(t, t0) = 1,5 (1 - exp(-(t - t0) / 100)), module constant
phi = lambda t, t0: 1.5 * (1 - math.exp(-(t - t0) / 100)) if t > t0 else 0.0
u = poteau(fonctions=(phi, lambda t0: 35e9))
print(f"{'utilisateur':17} [mm] : " + " ".join(f"{x * 1e3:8.3f}" for x in u))
print(f"contrôle COURBON à 10000 j : u0 (1 + phi) = {u0 * (1 + 2 * (1 - math.exp(-0.006 * 9972))) * 1e3:.3f} mm")
