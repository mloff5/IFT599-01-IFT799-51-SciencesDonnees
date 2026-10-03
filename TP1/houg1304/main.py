# -*- coding: utf-8 -*-
# =============================================================================
# IFT599 - Science de donnees
# TP1 : Comprehension et visualisation des donnees - Automne 2026
# Classe : IFT599-1
#
# Equipe :
#   Mael Hounkou                 (houg1304)
#   Mouhamadou Moustapha Gueye   (guem8961)
#   Cheikh Ahmadou Tidiane Ly    (lyzc3821)
#
# Fichier : main.py
# Role    : programme principal. Enchaine le chargement des donnees, la
#           Methode 1 (mesures quantitatives) et la Methode 2 (selection de
#           variables et visualisation), puis produit la synthese des resultats.
#
# Utilisation :
#   python main.py                        # execution complete
#   python main.py --etape m1             # Methode 1 seulement
#   python main.py --etape m2             # Methode 2 seulement
#   python main.py --n-variables 1000     # autre taille d'ensemble de depart
#   python main.py --donnees /chemin/vers/TCGA-PANCAN-HiSeq-801x20531
#
# Sources : voir l'en-tete de chaque module dans src/ et le fichier SOURCES.md.
# =============================================================================
"""Programme principal du TP1."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

RACINE = Path(__file__).resolve().parent
sys.path.insert(0, str(RACINE))

import config  # noqa: E402
from src import data_loader, methode1, methode2, synthese  # noqa: E402


def banniere():
    largeur = 74
    print("=" * largeur)
    print("IFT599 - Science de donnees - TP1".center(largeur))
    print("Comprehension et visualisation des donnees - Automne 2026".center(largeur))
    print("Classe : {}".format(config.CLASSE).center(largeur))
    print("-" * largeur)
    for nom, cip in config.EQUIPE:
        print("  {:<34s} ({})".format(nom, cip))
    print("=" * largeur)
    print()


def analyser_arguments():
    analyseur = argparse.ArgumentParser(
        description="TP1 IFT599 - separabilite et visualisation des profils genomiques.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    analyseur.add_argument(
        "--donnees", type=str, default=str(config.DOSSIER_DONNEES),
        help="dossier contenant data.csv et labels.csv")
    analyseur.add_argument(
        "--etape", choices=["tout", "m1", "m2"], default="tout",
        help="partie du TP a executer")
    analyseur.add_argument(
        "--n-variables", type=int, default=config.N_VARIABLES,
        help="taille de l'ensemble de donnees de depart")
    analyseur.add_argument(
        "--mode-selection", choices=["aleatoire", "variance"],
        default=config.MODE_SELECTION_INITIALE,
        help="mode de reduction dimensionnelle initiale")
    analyseur.add_argument(
        "--graine", type=int, default=config.GRAINE,
        help="graine aleatoire (reproductibilite)")
    analyseur.add_argument(
        "--sans-cache", action="store_true",
        help="recharge le fichier CSV au lieu d'utiliser le cache")
    analyseur.add_argument(
        "--sans-histogrammes", action="store_true",
        help="desactive les histogrammes marginaux dans les nuages de points")
    analyseur.add_argument(
        "--sortie", type=str, default=str(config.DOSSIER_RESULTATS),
        help="dossier de sortie des resultats")
    return analyseur.parse_args()


def main():
    arguments = analyser_arguments()
    banniere()

    # Application des options de ligne de commande a la configuration globale.
    config.N_VARIABLES = arguments.n_variables
    config.GRAINE = arguments.graine
    config.MODE_SELECTION_INITIALE = arguments.mode_selection
    if arguments.sans_histogrammes:
        config.AFFICHER_HISTOGRAMMES = False

    dossier_sortie = Path(arguments.sortie)
    config.DOSSIER_RESULTATS = dossier_sortie
    config.DOSSIER_TABLEAUX = dossier_sortie / "tableaux"
    config.DOSSIER_FIGURES = dossier_sortie / "figures"
    config.DOSSIER_TABLEAUX.mkdir(parents=True, exist_ok=True)
    config.DOSSIER_FIGURES.mkdir(parents=True, exist_ok=True)

    depart = time.time()

    # --- 1. Chargement et reduction dimensionnelle initiale ------------------
    print("[1/4] Chargement des donnees et construction de l'ensemble de depart")
    try:
        jeu = data_loader.charger(
            dossier_donnees=arguments.donnees,
            n_variables=arguments.n_variables,
            graine=arguments.graine,
            mode=arguments.mode_selection,
            utiliser_cache=not arguments.sans_cache,
        )
    except FileNotFoundError as erreur:
        print()
        print("ERREUR : " + str(erreur))
        return 1
    print(jeu.resume())
    chemin_variables = data_loader.ecrire_liste_variables(
        jeu, dossier_sortie / "variables_selectionnees.txt")
    print("  liste des variables retenues -> {}".format(chemin_variables.name))
    print()

    resultats_m1 = None
    resultats_m2 = None

    # --- 2. Methode 1 --------------------------------------------------------
    if arguments.etape in ("tout", "m1"):
        print("[2/4] Methode 1 - analyse de separabilite sans visualisation")
        resultats_m1 = methode1.executer(
            jeu,
            dossier_tableaux=config.DOSSIER_TABLEAUX,
            dossier_figures=config.DOSSIER_FIGURES,
        )
        nb_ok = int(resultats_m1["separabilite"]["bien_separees"].sum())
        total = len(resultats_m1["separabilite"])
        print("  {} combinaisons paire-de-classes bien separees sur {}"
              .format(nb_ok, total))
        print()

    # --- 3. Methode 2 --------------------------------------------------------
    if arguments.etape in ("tout", "m2"):
        print("[3/4] Methode 2 - selection de variables et visualisation")
        resultats_m2 = methode2.executer(
            jeu,
            dossier_tableaux=config.DOSSIER_TABLEAUX,
            dossier_figures=config.DOSSIER_FIGURES,
        )
        figures = sorted(config.DOSSIER_FIGURES.glob("*.png"))
        print("  {} figures produites dans {}".format(
            len(figures), config.DOSSIER_FIGURES))
        print()

    # --- 4. Synthese ---------------------------------------------------------
    print("[4/4] Redaction de la synthese des resultats")
    chemin_synthese = synthese.ecrire(
        jeu, resultats_m1, resultats_m2,
        chemin=dossier_sortie / "SYNTHESE.md")
    print("  -> {}".format(chemin_synthese))
    print()

    duree = time.time() - depart
    print("Execution terminee en {:.1f} s.".format(duree))
    print("Tableaux : {}".format(config.DOSSIER_TABLEAUX))
    print("Figures  : {}".format(config.DOSSIER_FIGURES))
    return 0


if __name__ == "__main__":
    sys.exit(main())
