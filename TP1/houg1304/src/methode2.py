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
# Fichier : src/methode2.py
# Role    : Methode 2 - analyse avec visualisation.
#           Selection automatique de paires de variables (originales, ACP) et
#           projection t-SNE, puis production des nuages de points.
#
# Sources : numpy, pandas, matplotlib ; scikit-learn pour t-SNE (cite dans
#           src/transformations.py). La selection de variables est une
#           implantation de l'equipe (voir src/selection.py).
# =============================================================================
"""Methode 2 : selection de paires de variables et visualisation.

Trois modalites de selection sont traitees, comme exige par l'enonce :

1. Variables originales   : selection automatique par la mesure Overlap()
                            (schema de la diapositive 37 du Theme 1).
                            Deux versions sont produites afin d'appuyer la
                            comparaison "transformation vs non-transformation" :
                            variables brutes et variables standardisees.
2. Variables transformees par ACP : 15 composantes principales generees, puis
                            selection de la meilleure paire avec le meme
                            algorithme.
3. t-SNE                  : la methode genere directement deux dimensions ;
                            aucune selection n'est necessaire. Un t-SNE global
                            fournit la projection valable pour toutes les
                            paires de classes, et un t-SNE restreint a chaque
                            couple de classes fournit la projection propre a
                            chaque paire.

Pour chaque modalite, on determine :
  * la meilleure paire de variables pour TOUTES les paires de classes
    (mesure d'evaluation : maximum des Overlap, critere minimax) ;
  * la meilleure paire de variables pour CHAQUE paire de classes
    (mesure d'evaluation : Overlap de la paire de classes consideree).

La distance utilisee pour toute la Methode 2 est la distance Euclidienne,
comme conseille dans l'enonce.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config  # noqa: E402

from .selection import (  # noqa: E402
    EvaluateurGlobal, EvaluateurPaireDeClasses, rechercher_meilleure_paire,
)
from .separabilite import details_paire  # noqa: E402
from .transformations import ACP, appliquer_transformation, tsne  # noqa: E402
from .visualisation import courbe_variance_acp, nuage_de_points  # noqa: E402


# -----------------------------------------------------------------------------
# Representations des donnees utilisees par la Methode 2
# -----------------------------------------------------------------------------
def construire_representations(jeu, verbeux=True):
    """Construit les espaces de variables dans lesquels se fait la selection."""
    representations = {}

    # --- 1a. Variables originales, sans transformation -----------------------
    representations["originales_brutes"] = {
        "matrice": appliquer_transformation(jeu.X, "brut"),
        "noms": list(jeu.noms_variables),
        "titre": "variables originales (brutes)",
        "prefixe": "originales_brutes",
    }

    # --- 1b. Variables originales standardisees ------------------------------
    representations["originales_zscore"] = {
        "matrice": appliquer_transformation(jeu.X, "zscore"),
        "noms": list(jeu.noms_variables),
        "titre": "variables originales (Z-score)",
        "prefixe": "originales_zscore",
    }

    # --- 2. Composantes principales -----------------------------------------
    if verbeux:
        print("  - ACP : generation de {} composantes principales..."
              .format(config.N_COMPOSANTES_ACP))
    acp = ACP(n_composantes=config.N_COMPOSANTES_ACP, centrer_reduire=True)
    scores_acp = acp.ajuster_transformer(jeu.X)
    representations["acp"] = {
        "matrice": scores_acp,
        "noms": acp.noms_composantes(),
        "titre": "composantes principales (ACP)",
        "prefixe": "acp",
        "modele": acp,
    }
    if verbeux:
        print("    variance expliquee par les {} composantes : {:.1f} %"
              .format(config.N_COMPOSANTES_ACP, acp.variance_cumulee_[-1] * 100.0))

    return representations


# -----------------------------------------------------------------------------
# Execution de la Methode 2
# -----------------------------------------------------------------------------
def executer(jeu, dossier_tableaux=None, dossier_figures=None, verbeux=True):
    dossier_tableaux = Path(dossier_tableaux or config.DOSSIER_TABLEAUX)
    dossier_figures = Path(dossier_figures or config.DOSSIER_FIGURES)
    dossier_tableaux.mkdir(parents=True, exist_ok=True)
    dossier_figures.mkdir(parents=True, exist_ok=True)

    classes = sorted(jeu.classes)
    paires_de_classes = jeu.paires_de_classes()
    lignes_resultats = []
    traces = []

    representations = construire_representations(jeu, verbeux=verbeux)

    # Figure d'appui : eboulis des valeurs propres de l'ACP.
    acp = representations["acp"]["modele"]
    courbe_variance_acp(acp.variance_expliquee_, acp.variance_cumulee_,
                        dossier_figures / "acp_variance_expliquee.png")

    # -------------------------------------------------------------------------
    # A. Modalites avec selection de variables (originales et ACP)
    # -------------------------------------------------------------------------
    for cle, representation in representations.items():
        M = representation["matrice"]
        noms = representation["noms"]
        titre_modalite = representation["titre"]
        prefixe = representation["prefixe"]

        if verbeux:
            print("  - selection sur les {} ({} variables)"
                  .format(titre_modalite, M.shape[1]))

        # --- A.1 Meilleure paire pour TOUTES les paires de classes -----------
        # Deux mesures d'evaluation agregees sont comparees (voir selection.py).
        for agregation in ("max", "moyenne"):
            evaluateur = EvaluateurGlobal(jeu.y, classes, agregation=agregation)
            resultat = rechercher_meilleure_paire(
                M, evaluateur,
                n_candidates=config.N_CANDIDATES,
                max_iterations=config.MAX_ITERATIONS,
                epsilon=config.EPSILON_ARRET,
            )
            j, k = resultat["variable_1"], resultat["variable_2"]
            nom_j, nom_k = noms[j], noms[k]
            traces.append(
                ("{} | toutes les paires de classes | agregation '{}'".format(
                    titre_modalite, agregation), resultat["trace"]))

            lignes_resultats.append({
                "modalite": titre_modalite,
                "cible": "toutes les paires de classes",
                "agregation": agregation,
                "variable_1": nom_j,
                "variable_2": nom_k,
                "mesure_evaluation": resultat["mesure"],
                "score": resultat["score"],
                "overlap_max": resultat["overlap_max"],
                "overlap_moyen": resultat["overlap_moyen"],
                "paires_bien_separees": "{}/{}".format(
                    resultat["nb_paires_bien_separees"], resultat["nb_paires"]),
                "nb_paires_evaluees": resultat["nb_paires_evaluees"],
                "nb_iterations": resultat["nb_iterations"],
                "strategie": resultat["strategie"],
            })

            nuage_de_points(
                M[:, [j, k]], jeu.y, classes,
                "Methode 2 - {} : meilleure paire pour les 5 classes".format(titre_modalite),
                nom_j, nom_k,
                dossier_figures / "{}_toutes_classes_{}.png".format(prefixe, agregation),
                sous_titre="Critere : {} = {:.3f}  |  Overlap moyen = {:.3f}  |  "
                           "{} paires de classes bien separees sur {}".format(
                               resultat["mesure"], resultat["score"],
                               resultat["overlap_moyen"],
                               resultat["nb_paires_bien_separees"],
                               resultat["nb_paires"]),
            )

        # --- A.2 Meilleure paire pour CHAQUE paire de classes ----------------
        for classe_1, classe_2 in paires_de_classes:
            evaluateur_paire = EvaluateurPaireDeClasses(jeu.y, classe_1, classe_2)
            resultat_paire = rechercher_meilleure_paire(
                M, evaluateur_paire,
                n_candidates=config.N_CANDIDATES,
                max_iterations=config.MAX_ITERATIONS,
                epsilon=config.EPSILON_ARRET,
            )
            j2, k2 = resultat_paire["variable_1"], resultat_paire["variable_2"]
            nom_j2, nom_k2 = noms[j2], noms[k2]
            traces.append(("{} | {} vs {}".format(titre_modalite, classe_1, classe_2),
                           resultat_paire["trace"]))

            lignes_resultats.append({
                "modalite": titre_modalite,
                "cible": "{} / {}".format(classe_1, classe_2),
                "agregation": "-",
                "variable_1": nom_j2,
                "variable_2": nom_k2,
                "mesure_evaluation": resultat_paire["mesure"],
                "score": resultat_paire["score"],
                "overlap_max": resultat_paire["score"],
                "overlap_moyen": resultat_paire["score"],
                "paires_bien_separees": "1/1" if resultat_paire["score"] < 1.0 else "0/1",
                "nb_paires_evaluees": resultat_paire["nb_paires_evaluees"],
                "nb_iterations": resultat_paire["nb_iterations"],
                "strategie": resultat_paire["strategie"],
            })

            masque = np.isin(jeu.y, [classe_1, classe_2])
            nuage_de_points(
                M[np.ix_(np.flatnonzero(masque), [j2, k2])],
                jeu.y[masque], [classe_1, classe_2],
                "Methode 2 - {} : {} vs {}".format(titre_modalite, classe_1, classe_2),
                nom_j2, nom_k2,
                dossier_figures / "{}_{}_{}.png".format(prefixe, classe_1, classe_2),
                sous_titre="Overlap = {:.3f} ({})".format(
                    resultat_paire["score"],
                    "classes bien separees" if resultat_paire["score"] < 1.0
                    else "superposition detectee"),
            )

    # -------------------------------------------------------------------------
    # B. t-SNE : les deux dimensions sont generees directement
    # -------------------------------------------------------------------------
    if verbeux:
        print("  - t-SNE global sur les {} patients...".format(jeu.n))

    X_standardise = appliquer_transformation(jeu.X, "zscore")
    projection_globale, perplexite_utilisee = tsne(
        X_standardise,
        n_composantes=2,
        perplexite=config.TSNE_PERPLEXITE,
        graine=config.GRAINE,
        n_iter=config.TSNE_N_ITER,
        learning_rate=config.TSNE_LEARNING_RATE,
    )

    overlaps_tsne = {}
    for classe_1, classe_2 in paires_de_classes:
        details = details_paire(projection_globale[jeu.y == classe_1],
                                projection_globale[jeu.y == classe_2])
        overlaps_tsne["{}-{}".format(classe_1, classe_2)] = details["overlap"]

    valeurs_tsne = np.array(list(overlaps_tsne.values()), dtype=np.float64)
    finies = valeurs_tsne[np.isfinite(valeurs_tsne)]
    nb_separees_tsne = int((valeurs_tsne < 1.0).sum())
    lignes_resultats.append({
        "modalite": "t-SNE (projection globale)",
        "cible": "toutes les paires de classes",
        "agregation": "max et moyenne",
        "variable_1": "t-SNE 1",
        "variable_2": "t-SNE 2",
        "mesure_evaluation": "aucune selection : les 2 dimensions sont generees",
        "score": float(valeurs_tsne.max()),
        "overlap_max": float(valeurs_tsne.max()),
        "overlap_moyen": float(finies.mean()) if finies.size else np.nan,
        "paires_bien_separees": "{}/{}".format(nb_separees_tsne, len(overlaps_tsne)),
        "nb_paires_evaluees": 0,
        "nb_iterations": 0,
        "strategie": "generation directe de 2 dimensions (perplexite {:.0f})".format(
            perplexite_utilisee),
    })

    nuage_de_points(
        projection_globale, jeu.y, classes,
        "Methode 2 - t-SNE : projection globale des 5 classes",
        "t-SNE 1", "t-SNE 2",
        dossier_figures / "tsne_toutes_classes.png",
        sous_titre="max des Overlap = {:.3f}  |  Overlap moyen = {:.3f}  |  "
                   "{} paires de classes bien separees sur {}".format(
                       valeurs_tsne.max(),
                       finies.mean() if finies.size else np.nan,
                       nb_separees_tsne, len(overlaps_tsne)),
    )

    # t-SNE restreint a chaque paire de classes.
    for classe_1, classe_2 in paires_de_classes:
        if verbeux:
            print("    t-SNE {} vs {}".format(classe_1, classe_2))
        masque = np.isin(jeu.y, [classe_1, classe_2])
        projection, perplexite_paire = tsne(
            X_standardise[masque],
            n_composantes=2,
            perplexite=config.TSNE_PERPLEXITE,
            graine=config.GRAINE,
            n_iter=config.TSNE_N_ITER,
            learning_rate=config.TSNE_LEARNING_RATE,
        )
        y_paire = jeu.y[masque]
        details = details_paire(projection[y_paire == classe_1],
                                projection[y_paire == classe_2])

        lignes_resultats.append({
            "modalite": "t-SNE (projection par paire)",
            "cible": "{} / {}".format(classe_1, classe_2),
            "agregation": "-",
            "variable_1": "t-SNE 1",
            "variable_2": "t-SNE 2",
            "mesure_evaluation": "Overlap({}, {})".format(classe_1, classe_2),
            "score": details["overlap"],
            "overlap_max": details["overlap"],
            "overlap_moyen": details["overlap"],
            "paires_bien_separees": "1/1" if details["overlap"] < 1.0 else "0/1",
            "nb_paires_evaluees": 0,
            "nb_iterations": 0,
            "strategie": "generation directe de 2 dimensions (perplexite {:.0f})".format(
                perplexite_paire),
        })

        nuage_de_points(
            projection, y_paire, [classe_1, classe_2],
            "Methode 2 - t-SNE : {} vs {}".format(classe_1, classe_2),
            "t-SNE 1", "t-SNE 2",
            dossier_figures / "tsne_{}_{}.png".format(classe_1, classe_2),
            sous_titre="Overlap = {:.3f} ({})".format(
                details["overlap"],
                "classes bien separees" if details["overlap"] < 1.0
                else "superposition detectee"),
        )

    tableau = pd.DataFrame(lignes_resultats)
    tableau.to_csv(dossier_tableaux / "methode2_selection.csv", index=False)

    # Trace detaillee du deroulement de la selection automatique (schema p.37).
    chemin_trace = dossier_tableaux / "methode2_trace_selection.txt"
    with open(chemin_trace, "w", encoding="utf-8") as f:
        f.write(config.entete_texte(
            "Methode 2 - trace de la selection automatique de variables "
            "(schema de la diapositive 37 du Theme 1)"))
        for titre, lignes in traces:
            f.write("\n### {}\n".format(titre))
            for ligne in lignes:
                f.write("  " + ligne + "\n")

    return {
        "selection": tableau,
        "overlaps_tsne": overlaps_tsne,
        "acp": acp,
    }
