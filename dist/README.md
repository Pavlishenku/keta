# keta 1.0.0 — installation

## Paquet Python

```sh
pip install keta
python -c "import keta; print(keta.__version__)"
```

Des paquets précompilés sont fournis pour Python 3.10 à 3.13, sous :

- Linux x86-64, glibc 2.28 ou plus (Debian 10, Ubuntu 20.04, RHEL 8 et suivants) ;
- Windows 64 bits.

Le paquet n'a aucune dépendance Python. Il installe aussi le programme `keta` en ligne de commande.

## Compilation depuis les sources

```sh
pip install .
```

À lancer à la racine des sources. Il faut Python 3.10 ou plus et un compilateur C++20 : GCC 11 ou
plus sous Linux, Visual Studio 2022 (charge de travail « Développement Desktop en C++ ») sous
Windows. CMake et Ninja sont installés par `pip`. La compilation dure de 10 à 20 minutes.

macOS n'est pas pris en charge : sa bibliothèque standard ne fournit pas encore la lecture des
nombres réels par `std::from_chars`, que le noyau utilise.

## Interface graphique (Windows)

`keta_graphique-0.9.0-windows-x64.zip`, dans ce répertoire, contient keta graphique, l'interface
graphique de keta, pour Windows 10 ou 11 (64 bits). L'application est portable : décompresser
l'archive dans un dossier quelconque, puis lancer `keta_graphique.bat` (ou `bin\keta_graphique.exe`).
Les bibliothèques Qt qu'elle utilise et son manuel (`keta_graphique_manuel.pdf`) sont dans
l'archive. Le noyau de calcul est compris dans l'application : Python n'est nécessaire que pour y
exécuter des scripts.

L'application n'est pas signée : Windows SmartScreen peut afficher « Windows a protégé votre
ordinateur » ; cliquer sur « Informations complémentaires », puis « Exécuter quand même ». Sur une
machine sans pilote OpenGL 3.3 (bureau à distance, machine virtuelle), définir la variable
d'environnement `QT_OPENGL=software`.

## Ligne de commande

```sh
keta                                        # liste des commandes
keta run modele.json --out resultats.json   # statique d'un modèle enregistré par Model.to_json()
```

## Serveur de documentation

```sh
pip install "keta[mcp]"
keta-mcp
```

Le serveur donne accès, par le protocole MCP, au texte des deux manuels, à leurs exemples et à
l'API du module.
