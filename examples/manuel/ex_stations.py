# Sections d'étude : découpage des barres, lecture des efforts, contraintes et
# « contrainte généralisée » (combinaison linéaire d'efforts) calculée en Python.
import keta

def z(x):                                         # supprime le zéro négatif d'affichage
    return 0.0 if abs(x) < 1e-6 else x

L, q = 12.0, -40e3
m = keta.Model()
a, b = m.add_node(0, 0, 0), m.add_node(L, 0, 0)
poutre = m.add_beam(a, b, m.add_section_rect(0.4, 1.0), m.concrete())
m.add_support(a, dx=True, dy=True, dz=True, rx=True)
m.add_support(b, dy=True, dz=True)
m.add_udl(poutre, wz=q)
m.add_point_load(b, fx=-500e3)                    # compression de 500 kN

cas = m.analyze(stations=7).case(1)               # 7 sections d'étude par barre (pas L/6)
A, I, v = 0.4 * 1.0, 0.4 * 1.0**3 / 12, 0.5
print("  x/L      x [m]      N [kN]    Vz [kN]    My [kN.m]  sig_inf [MPa]  sig_sup [MPa]")
for x, N, Vy, Vz, T, My, Mz in cas.bar_forces(poutre):
    s_inf = N / A + My * (-v) / I                 # sigma = N/A + My z / Iy
    s_sup = N / A + My * (+v) / I
    print(f"{x:6.4f} {x * L:9.3f} {N / 1e3:11.2f} {z(Vz) / 1e3:10.2f} {z(My) / 1e3:12.2f} "
          f"{s_inf / 1e6:14.3f} {s_sup / 1e6:14.3f}")
# contrainte généralisée : C = 0,25 N + 0,33 My + 10 (combinaison quelconque d'efforts)
C = [0.25 * s[1] + 0.33 * s[5] + 10.0 for s in cas.bar_forces(poutre)]
print("contrainte généralisée C(x) :", " ".join(f"{c / 1e3:.2f}" for c in C), "(x 1000)")
print("efforts d'extrémité (repère local, 12 composantes) :")
print(" ", " ".join(f"{f / 1e3:.1f}" for f in cas.end_forces(poutre)))
