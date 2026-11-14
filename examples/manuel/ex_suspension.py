# Pont suspendu : un tablier de 21 travées de 10 m pendu à un câble porteur entre deux selles
# distantes de 240 m. Recherche de forme du câble sous les charges des suspentes, assemblage, puis
# suspentes à leur longueur à vide, correction en deux étapes des longueurs à vide du câble (le tablier,
# raide, reprend une part de la charge).
import keta

def pont():
    m = keta.Model()
    acier = m.add_material(200e9, 0.3, rho=0.0)
    fil = m.add_material(195e9, 0.3, rho=7850.0)
    tige = m.add_material(160e9, 0.3, rho=0.0)
    s_tablier = m.add_section_props(0.5, 0.8, 5.0, 1.0)
    s_cable = m.add_section_props(0.30, 1e-6, 1e-6, 1e-6)
    s_suspente = m.add_section_props(0.01, 1e-8, 1e-8, 1e-8)
    nh, x0, dx, q = 21, 20.0, 10.0, 150e3
    noeuds = [m.add_node(x0 - dx / 2, 0, 5.0)] + [m.add_node(x0 + dx * i, 0, 5.0) for i in range(nh)] \
        + [m.add_node(x0 + dx * (nh - 1) + dx / 2, 0, 5.0)]
    barres = [m.add_beam(a, b, s_tablier, acier) for a, b in zip(noeuds, noeuds[1:])]
    m.add_support(noeuds[0], dx=True, dy=True, rx=True)      # le tablier est pendu, pas appuyé
    m.add_support(noeuds[-1], dy=True)
    for b in barres:
        m.add_udl(b, wz=-q, case=1)
    donnees = dict(left=(0, 0, 40), right=(240, 0, 40), hanger_x=[x0 + dx * i for i in range(nh)],
                   hanger_load=[q * dx] * nh, hanger_deck_nodes=noeuds[1:-1], cable_material=fil,
                   cable_section=s_cable, hanger_material=tige, hanger_section=s_suspente, target_sag=25.0)
    return m, donnees, q * dx * nh

# 1. recherche de forme seule : tension horizontale et longueur à vide du câble
m, donnees, poids = pont()
forme = keta.form_find_cable(donnees["left"], donnees["right"], donnees["hanger_x"], donnees["hanger_load"],
                             ea=195e9 * 0.30, w=7850.0 * 0.30 * 9.80665, target_sag=25.0)
print(f"forme : H = {forme['H'] / 1e6:.2f} MN, flèche {forme['sag']:.3f} m, longueur à vide {forme['total_l0']:.3f} m")

# 2. assemblage dans le modèle et équilibre sous le poids du tablier
sm = m.build_suspension(**donnees)
r = m.cable_static(case=1)
ecart = max(abs(r.displacement(n)[2]) for n in sm["cable_nodes"][1:-1])
print(f"assemblé : {len(sm['cables'])} tronçons de câble, {len(sm['hanger_bars'])} suspentes ; "
      f"réaction verticale {r.sum_reactions()[2] / 1e6:.3f} MN pour {poids / 1e6:.3f} MN de tablier ; "
      f"câble écarté de son profil de {ecart * 1e3:.2f} mm")

# 3. les mêmes suspentes à leur longueur à vide : leur défaut de longueur va dans le cas permanent
m1, donnees1, _ = pont()
sm1 = m1.build_suspension(**donnees1, hanger_fit_case=1)
r1 = m1.cable_static(case=1)
fleche = min(r.displacement(n)[2] for n in donnees["hanger_deck_nodes"])
fleche1 = min(r1.displacement(n)[2] for n in donnees1["hanger_deck_nodes"])
print(f"suspentes : posées sans contrainte, le tablier descend de {-fleche * 1e3:.1f} mm ; à leur longueur à vide "
      f"(de {min(sm1['hanger_l0']):.3f} à {max(sm1['hanger_l0']):.3f} m), de {-fleche1 * 1e3:.1f} mm")

# 4. correction en deux étapes : les longueurs à vide ramènent le câble sur son profil
m2, donnees2, _ = pont()
u = m2.solve_uelm(1, **donnees2, tol=1e-10)
r2 = m2.cable_static(case=1)
ecart2 = max(abs(r2.displacement(n)[2]) for n in u["cable_nodes"][1:-1])
print(f"deux étapes : écart {u['deviation_stage1'] * 1e3:.2f} mm avant correction, {ecart2 * 1e3:.4f} mm après, "
      f"{u['iterations']} itérations ; longueur à vide totale {sum(u['l0']):.3f} m")
assert ecart2 < 1e-7 and abs(r.sum_reactions()[2] - poids - 7850.0 * 0.30 * 9.80665 * forme["total_l0"]) < 1e-6 * poids
