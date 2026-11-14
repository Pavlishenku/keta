#!/usr/bin/env python3
"""Exécute tous les exemples du manuel d'utilisation de l'API Python de keta
et vérifie leurs sorties.

Chaque exemple ``*.py`` de ce répertoire est un script autonome : il construit
un modèle, lance un calcul, imprime ses résultats et contrôle lui-même les
valeurs qui ont une forme fermée (``assert``). Ce script :

  * exécute chaque exemple dans un processus séparé, le module ``keta`` étant
    pris dans le répertoire donné par ``--module`` (par défaut la variable
    d'environnement ``KETA_PYTHONPATH``, sinon le module installé, sinon
    ``<dépôt>/build/python``) ;
  * exige un code de retour nul (aucune assertion en échec, aucune exception) ;
  * compare la sortie standard, caractère pour caractère, à la sortie de
    référence ``sorties/<exemple>.txt`` — celle que le manuel reproduit.

Usage ::

    python3 executer_tous.py                                # vérification
    python3 executer_tous.py --maj                          # (ré)écrit les références
    python3 executer_tous.py --module build/python ex_console
    python3 executer_tous.py --dossier ../exemples_en       # les exemples d'un autre répertoire

Code de sortie : 0 si tous les exemples passent, 1 sinon.
"""

from __future__ import annotations

import argparse
import difflib
import os
import subprocess
import sys
import time
from pathlib import Path

ICI = Path(__file__).resolve().parent
DEPOT = ICI.parents[2]


def exemples(dossier: Path, filtres: list[str]) -> list[Path]:
    tous = sorted(p for p in dossier.glob("*.py") if p.name != Path(__file__).name)
    if filtres:
        tous = [p for p in tous if any(f in p.stem for f in filtres)]
    return tous


def executer(script: Path, module: Path, delai: float) -> tuple[int, str, str, float]:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(module) + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    env["PYTHONHASHSEED"] = "0"
    env["PYTHONIOENCODING"] = "utf-8"
    t0 = time.perf_counter()
    p = subprocess.run([sys.executable, str(script)], cwd=str(script.parent), env=env, capture_output=True,
                       text=True, timeout=delai)
    return p.returncode, p.stdout, p.stderr, time.perf_counter() - t0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("filtres", nargs="*", help="sous-chaînes des noms d'exemples à exécuter (tous par défaut)")
    ap.add_argument("--module", default=os.environ.get("KETA_PYTHONPATH"),
                    help="répertoire contenant le paquet keta (défaut : le module installé, sinon build/python)")
    ap.add_argument("--maj", action="store_true", help="écrire les sorties obtenues comme références")
    ap.add_argument("--dossier", default=str(ICI),
                    help="dossier des exemples (défaut : celui de ce script ; exemples_en pour l'anglais)")
    ap.add_argument("--delai", type=float, default=600.0, help="délai maximal par exemple [s]")
    a = ap.parse_args()
    if a.module is None:
        try:
            import keta
            a.module = str(Path(keta.__file__).resolve().parents[1])
        except ImportError:
            a.module = str(DEPOT / "build" / "python")
    module = Path(a.module).resolve()
    if not (module / "keta" / "__init__.py").exists():
        print(f"paquet keta introuvable dans {module}", file=sys.stderr)
        return 1
    dossier = Path(a.dossier).resolve()
    SORTIES = dossier / "sorties"
    SORTIES.mkdir(exist_ok=True)
    liste = exemples(dossier, a.filtres)
    echecs = 0
    for s in liste:
        try:
            code, out, err, dt = executer(s, module, a.delai)
        except subprocess.TimeoutExpired:
            print(f"DÉLAI   {s.name}")
            echecs += 1
            continue
        ref = SORTIES / (s.stem + ".txt")
        if code != 0:
            print(f"ÉCHEC   {s.name}  (code {code}, {dt:.1f} s)")
            print("        " + "\n        ".join((err or out).strip().splitlines()[-12:]))
            echecs += 1
            continue
        if a.maj:
            ref.write_text(out, encoding="utf-8")
            print(f"ÉCRIT   {s.name}  ({dt:.1f} s)")
            continue
        if not ref.exists():
            print(f"SANS RÉFÉRENCE  {s.name}  (lancer avec --maj)")
            echecs += 1
            continue
        attendu = ref.read_text(encoding="utf-8")
        if out != attendu:
            print(f"ÉCART   {s.name}")
            for ligne in list(difflib.unified_diff(attendu.splitlines(), out.splitlines(), "référence",
                                                   "obtenu", lineterm=""))[:30]:
                print("        " + ligne)
            echecs += 1
        else:
            print(f"OK      {s.name}  ({dt:.1f} s)")
    print(f"{len(liste) - echecs}/{len(liste)} exemple(s) conforme(s)")
    return 1 if echecs else 0


if __name__ == "__main__":
    sys.exit(main())
