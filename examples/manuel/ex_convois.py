# Charges d'exploitation routières sur une ligne d'influence : modèle LM1 par voie,
# systèmes A(l), Bc et Mc120 et répartition transversale de Courbon.
import keta

# Travée indépendante de 30 m : moment à mi-portée (ligne triangulaire exacte, échantillonnée)
L = 30.0
m = keta.Model()
n = [m.add_node(L * i / 12, 0, 0) for i in range(13)]
b = [m.add_beam(n[i], n[i + 1], m.add_section_props(5.0, 3.0, 30.0, 4.0), m.concrete()) for i in range(12)]
m.add_support(n[0], dx=True, dy=True, dz=True, rx=True)
m.add_support(n[-1], dy=True, dz=True)
lig = m.influence_line(b, dict(kind="bar_force", bar=b[5], x=2.5, comp=4), per_bar=21)
print(f"ordonnée maximale |eta| = {abs(lig.point_extreme()[1]):.4f} m (L/4 = {L / 4:.4f} m)")

# EN 1991-2, LM1 : voies de 3 m sur une chaussée de 9 m
nv, larg, reste = keta.en_lane_division(9.0)
print(f"chaussée de 9 m : {nv} voies de {larg} m, aire résiduelle {reste} m")
for rang in (1, 2, 3):
    res = keta.lm1_lane_extreme(lig, rang, larg)
    print(f"voie {rang} : tandem {keta.lm1_tandem(rang)['axles'][0][1] / 1e3:.0f} kN/essieu, "
          f"UDL {keta.lm1_udl(rang, larg) / 1e3:.1f} kN/m -> M = {min(res[:2]) / 1e3:9.1f} kN.m")

# Fascicule 61 titre II : A(l) sur la longueur chargée la plus défavorable, Bc, Mc120
dens = lambda Lc: keta.fasc61_a_lane_load(Lc, 1, 3.5) * 2      # 2 voies, pont de 1re classe [N/m]
val, Lc, zones = keta.intensity_udl_extreme(lig, dens, sign=-1)
print(f"A(l) : longueur chargée {Lc:.2f} m, A(30) = {keta.fasc61_a_of_l(30.0):.1f} N/m2 -> M = {val / 1e3:.1f} kN.m")
G, S = 30 * 150e3, 2 * 60 * 9.80665e3                            # poids propre et charge Bc maximale
delta = keta.fasc61_delta(L, G, S)
mx, mn, rmx, rmn = keta.convoy_extreme(lig, keta.fasc61_bc(), 3000)
print(f"Bc (une file) : M = {mn / 1e3:.1f} kN.m, tête du convoi à s = {rmn:.2f} m ; "
      f"majoration dynamique delta = {delta:.4f}")
mx, mn, rmx, rmn = keta.convoy_extreme(lig, keta.fasc61_mc120(), 3000)
print(f"Mc120 : M = {mn / 1e3:.1f} kN.m")

# Répartition transversale de Courbon : 4 poutres, charge unité à y = 3 m
K = keta.courbon_coefficients([-4.5, -1.5, 1.5, 4.5], y=3.0)
print("coefficients de Courbon :", [round(k, 4) for k in K])
