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
# Fichier : src/methode1.py
# Role    : Methode 1 - analyse de separabilite sans visualisation.
#           Les mesures sont calculees sur TOUTES les variables de l'ensemble
#           de donnees de depart et presentees sous forme de tableaux.
#
# Sources : numpy, pandas. Formules issues de l'enonce du TP1.
# =============================================================================
"""Methode 1 : mesures quantitatives de separation presentees en tableaux.

Pour chaque combinaison (transformation des donnees x metrique de distance),
on calcule :
  * la distance intra-classe de chacune des 5 classes ;
  * pour chacune des 10 paires mutuelles de classes : la distance inter-classe
    et l'indicateur de superposition Overlap().

Les resultats sont enregistres en CSV (pour le rapport) et resumes en tableaux
Markdown directement reutilisables.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config  # noqa: E402

from .distances import preparer_espace  # noqa: E402
from .separabilite import toutes_les_cohesions, toutes_les_paires  # noqa: E402
from .transformations import (  # noqa: E402
    ETIQUETTES_TRANSFORMATIONS, appliquer_transformation,
)
from .visualisation import carte_de_chaleur_overlap  # noqa: E402


def _matrice_overlap(resultats_paires, classes):
    """Construit la matrice carree symetrique des valeurs de Overlap."""
    classes = sorted(classes)
    index = {c: i for i, c in enumerate(classes)}
    matrice = np.full((len(classes), len(classes)), np.nan)
    for ligne in resultats_paires:
        i, j = index[ligne["classe_1"]], index[ligne["classe_2"]]
        matrice[i, j] = ligne["overlap"]
        matrice[j, i] = ligne["overlap"]
    return matrice, classes


def executer(jeu, transformations=None, metriques=None,
             dossier_tableaux=None, dossier_figures=None, verbeux=True):
    """Execute la Methode 1 pour toutes les combinaisons demandees.

    Retour
    ------
    dictionnaire contenant les tableaux (DataFrame) et les diagnostics.
    """
    transformations = transformations or config.TRANSFORMATIONS_METHODE1
    metriques = metriques or config.METRIQUES_METHODE1
    dossier_tableaux = Path(dossier_tableaux or config.DOSSIER_TABLEAUX)
    dossier_figures = Path(dossier_figures or config.DOSSIER_FIGURES)
    dossier_tableaux.mkdir(parents=True, exist_ok=True)
    dossier_figures.mkdir(parents=True, exist_ok=True)

    classes = sorted(jeu.classes)
    lignes_cohesion = []
    lignes_paires = []
    diagnostics = []

    for nom_transformation in transformations:
        X = appliquer_transformation(jeu.X, nom_transformation)

        for nom_metrique in metriques:
            if verbeux:
                print("  - transformation '{}' + distance {}"
                      .format(nom_transformation, nom_metrique))

            Z, infos = preparer_espace(
                X, nom_metrique,
                methode_inversion=config.METHODE_INVERSION_COVARIANCE,
                lambda_ridge=config.LAMBDA_RIDGE,
            )
            if nom_metrique == "mahalanobis":
                diagnostics.append({
                    "transformation": nom_transformation,
                    "metrique": nom_metrique,
                    **{k: v for k, v in infos.items()},
                })

            # --- Tableau 1 : cohesion (distance intra-classe) ----------------
            cohesions = toutes_les_cohesions(Z, jeu.y, classes)
            for classe, valeur in cohesions.items():
                lignes_cohesion.append({
                    "transformation": nom_transformation,
                    "metrique": nom_metrique,
                    "classe": classe,
                    "effectif": int(np.sum(jeu.y == classe)),
                    "dist_intra": valeur,
                })

            # --- Tableau 2 : separation et Overlap par paire -----------------
            resultats = toutes_les_paires(Z, jeu.y, classes)
            for ligne in resultats:
                lignes_paires.append({
                    "transformation": nom_transformation,
                    "metrique": nom_metrique,
                    "classe_1": ligne["classe_1"],
                    "classe_2": ligne["classe_2"],
                    "dist_intra_C1": ligne["dist_intra_C1"],
                    "dist_intra_C2": ligne["dist_intra_C2"],
                    "dist_C1_vers_C2": ligne["dist_C1_vers_C2"],
                    "dist_C2_vers_C1": ligne["dist_C2_vers_C1"],
                    "dist_inter": ligne["dist_inter"],
                    "overlap": ligne["overlap"],
                    "bien_separees": ligne["bien_separees"],
                })

            # --- Figure d'appui : carte de chaleur des Overlap ---------------
            matrice, ordre = _matrice_overlap(resultats, classes)
            carte_de_chaleur_overlap(
                matrice, ordre,
                "Methode 1 - Overlap ({}, distance {})".format(
                    ETIQUETTES_TRANSFORMATIONS.get(nom_transformation, nom_transformation),
                    nom_metrique),
                dossier_figures / "methode1_overlap_{}_{}.png".format(
                    nom_transformation, nom_metrique),
            )

    tableau_cohesion = pd.DataFrame(lignes_cohesion)
    tableau_paires = pd.DataFrame(lignes_paires)
    tableau_diagnostics = pd.DataFrame(diagnostics)

    tableau_cohesion.to_csv(dossier_tableaux / "methode1_cohesion.csv", index=False)
    tableau_paires.to_csv(dossier_tableaux / "methode1_separabilite.csv", index=False)
    if not tableau_diagnostics.empty:
        tableau_diagnostics.to_csv(
            dossier_tableaux / "methode1_diagnostic_covariance.csv", index=False)

    # Tableau synthetique : une ligne par paire de classes, une colonne par
    # combinaison transformation x metrique. C'est le tableau le plus lisible
    # pour le rapport.
    pivot = tableau_paires.copy()
    pivot["paire"] = pivot["classe_1"] + " / " + pivot["classe_2"]
    pivot["combinaison"] = pivot["transformation"] + " + " + pivot["metrique"]
    tableau_pivot = pivot.pivot(index="paire", columns="combinaison", values="overlap")
    tableau_pivot.to_csv(dossier_tableaux / "methode1_overlap_synthese.csv")

    return {
        "cohesion": tableau_cohesion,
        "separabilite": tableau_paires,
        "synthese_overlap": tableau_pivot,
        "diagnostics": tableau_diagnostics,
    }
