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
# Fichier : src/synthese.py
# Role    : mise en forme des resultats (tableaux Markdown) destines a etre
#           repris directement dans le rapport.
#
# Sources : pandas.
# =============================================================================
"""Generation du fichier de synthese des resultats."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config  # noqa: E402


def _tableau_markdown(df, flottants=4):
    """Convertit un DataFrame en tableau Markdown."""
    if df is None or len(df) == 0:
        return "_(aucun resultat)_\n"
    copie = df.copy()
    for colonne in copie.columns:
        if pd.api.types.is_float_dtype(copie[colonne]):
            copie[colonne] = copie[colonne].map(
                lambda v: "inf" if not np.isfinite(v) else "{:.{p}f}".format(v, p=flottants))
    entetes = [str(c) for c in copie.columns]
    lignes = ["| " + " | ".join(entetes) + " |",
              "|" + "|".join(["---"] * len(entetes)) + "|"]
    for _, ligne in copie.iterrows():
        lignes.append("| " + " | ".join(str(v) for v in ligne.to_list()) + " |")
    return "\n".join(lignes) + "\n"


def _entete_markdown():
    membres = "\n".join("- {} ({})".format(nom, cip) for nom, cip in config.EQUIPE)
    return (
        "# IFT599 - Science de donnees - TP1\n"
        "## Comprehension et visualisation des donnees (Automne 2026)\n\n"
        "**Classe :** {classe}\n\n"
        "**Equipe :**\n\n{membres}\n\n"
        "---\n\n"
        "> Ce fichier rassemble tous les resultats numeriques produits par le\n"
        "> programme. Les tableaux sont directement reutilisables dans le\n"
        "> rapport. Les figures correspondantes se trouvent dans\n"
        "> `resultats/figures/`.\n\n"
    ).format(classe=config.CLASSE, membres=membres)


def ecrire(jeu, resultats_m1, resultats_m2, chemin=None):
    """Ecrit le fichier resultats/SYNTHESE.md."""
    chemin = Path(chemin or (config.DOSSIER_RESULTATS / "SYNTHESE.md"))
    chemin.parent.mkdir(parents=True, exist_ok=True)

    morceaux = [_entete_markdown()]

    # --- 1. Ensemble de donnees ---------------------------------------------
    morceaux.append("## 1. Ensemble de donnees de depart\n\n")
    morceaux.append("```\n" + jeu.resume() + "\n```\n\n")
    morceaux.append(
        "La reduction de dimensionnalite initiale a consiste a retirer les "
        "variables constantes puis a tirer au hasard (graine {}) {} variables "
        "parmi les {} variables non constantes. La liste complete des variables "
        "retenues se trouve dans `resultats/variables_selectionnees.txt`.\n\n"
        .format(config.GRAINE, jeu.p, jeu.p_non_constant))

    # --- 2. Methode 1 --------------------------------------------------------
    if resultats_m1 is not None:
        morceaux.append("## 2. Methode 1 - separabilite sans visualisation\n\n")
        morceaux.append(
            "Mesures calculees sur **toutes** les {} variables de l'ensemble de "
            "depart.\n\n".format(jeu.p))

        morceaux.append("### 2.1 Synthese des valeurs de Overlap\n\n")
        morceaux.append(
            "Une valeur inferieure a 1 indique que les deux classes sont bien "
            "separees.\n\n")
        synthese = resultats_m1["synthese_overlap"].reset_index()
        morceaux.append(_tableau_markdown(synthese))
        morceaux.append("\n")

        morceaux.append("### 2.2 Cohesion (distance intra-classe)\n\n")
        cohesion = resultats_m1["cohesion"].pivot(
            index="classe", columns=["transformation", "metrique"], values="dist_intra")
        cohesion.columns = ["{} + {}".format(a, b) for a, b in cohesion.columns]
        morceaux.append(_tableau_markdown(cohesion.reset_index()))
        morceaux.append("\n")

        morceaux.append("### 2.3 Detail par combinaison\n\n")
        for (transformation, metrique), groupe in resultats_m1["separabilite"].groupby(
                ["transformation", "metrique"], sort=False):
            morceaux.append("#### {} + distance {}\n\n".format(transformation, metrique))
            colonnes = ["classe_1", "classe_2", "dist_intra_C1", "dist_intra_C2",
                        "dist_inter", "overlap", "bien_separees"]
            morceaux.append(_tableau_markdown(groupe[colonnes].reset_index(drop=True)))
            nb_ok = int(groupe["bien_separees"].sum())
            morceaux.append("\nPaires bien separees : **{}/{}**\n\n".format(
                nb_ok, len(groupe)))

        if not resultats_m1["diagnostics"].empty:
            morceaux.append("### 2.4 Diagnostic de la matrice de covariance "
                            "(distance de Mahalanobis)\n\n")
            morceaux.append(_tableau_markdown(
                resultats_m1["diagnostics"].reset_index(drop=True), flottants=3))
            morceaux.append(
                "\nAvec n = {} patients et p = {} variables, la covariance "
                "empirique est de rang au plus {} : elle est donc singuliere. "
                "La regularisation retenue est indiquee ci-dessus.\n\n"
                "Remarque utile pour la discussion : la distance de Mahalanobis "
                "calculee avec la covariance exacte serait rigoureusement "
                "invariante a une transformation affine par variable (le "
                "Z-score et le min-max en font partie). Les ecarts observes "
                "entre les colonnes `brut`, `zscore` et `minmax` en Mahalanobis "
                "proviennent donc uniquement de la regularisation, qui n'est "
                "pas invariante a l'echelle des variables.\n\n"
                .format(jeu.n, jeu.p, min(jeu.n - 1, jeu.p)))

    # --- 3. Methode 2 --------------------------------------------------------
    if resultats_m2 is not None:
        morceaux.append("## 3. Methode 2 - selection de variables et visualisation\n\n")
        tableau = resultats_m2["selection"]

        morceaux.append("### 3.1 Meilleure paire de variables pour "
                        "**toutes** les paires de classes\n\n")
        morceaux.append(
            "Deux mesures d'evaluation agregees sont comparees : le maximum des "
            "Overlap (critere minimax, le plus exigeant) et leur moyenne "
            "(qualite d'ensemble de l'affichage).\n\n")
        globales = tableau[tableau["cible"] == "toutes les paires de classes"]
        colonnes = ["modalite", "agregation", "variable_1", "variable_2",
                    "overlap_max", "overlap_moyen", "paires_bien_separees",
                    "nb_paires_evaluees", "nb_iterations"]
        morceaux.append(_tableau_markdown(globales[colonnes].reset_index(drop=True)))
        morceaux.append("\n")

        morceaux.append("### 3.2 Meilleure paire de variables pour "
                        "**chaque** paire de classes\n\n")
        for modalite, groupe in tableau[
                tableau["cible"] != "toutes les paires de classes"].groupby(
                    "modalite", sort=False):
            morceaux.append("#### {}\n\n".format(modalite))
            colonnes = ["cible", "variable_1", "variable_2", "score",
                        "nb_paires_evaluees", "nb_iterations"]
            sous_tableau = groupe[colonnes].rename(
                columns={"cible": "paire de classes", "score": "Overlap"})
            morceaux.append(_tableau_markdown(sous_tableau.reset_index(drop=True)))
            nb_ok = int((groupe["score"] < 1.0).sum())
            morceaux.append("\nPaires bien separees : **{}/{}**\n\n".format(
                nb_ok, len(groupe)))

        acp = resultats_m2.get("acp")
        if acp is not None:
            morceaux.append("### 3.3 ACP - variance expliquee\n\n")
            tableau_acp = pd.DataFrame({
                "composante": acp.noms_composantes(),
                "valeur_propre": acp.valeurs_propres_,
                "variance_expliquee_%": acp.variance_expliquee_ * 100.0,
                "variance_cumulee_%": acp.variance_cumulee_ * 100.0,
            })
            morceaux.append(_tableau_markdown(tableau_acp, flottants=3))
            morceaux.append("\n")

    # --- 4. Pistes d'analyse -------------------------------------------------
    morceaux.append("## 4. Elements de comparaison attendus dans le rapport\n\n")
    morceaux.append(
        "Les tableaux ci-dessus permettent de traiter les trois comparaisons "
        "exigees au point 3 de la section 2.2 de l'enonce :\n\n"
        "1. **Transformation vs non-transformation** : comparer les colonnes "
        "`brut`, `zscore` et `minmax` du tableau 2.1, ainsi que les modalites "
        "`variables originales (brutes)` et `variables originales (Z-score)` "
        "de la section 3.\n"
        "2. **Euclidienne vs Mahalanobis** : comparer, a transformation egale, "
        "les colonnes `+ euclidienne` et `+ mahalanobis` du tableau 2.1, en "
        "tenant compte du diagnostic de la section 2.4.\n"
        "3. **ACP vs t-SNE** : comparer les lignes correspondantes des tableaux "
        "3.1 et 3.2 et les figures `acp_*.png` et `tsne_*.png`.\n\n")
    morceaux.append("---\n\nFichier genere automatiquement par `main.py`.\n")

    with open(chemin, "w", encoding="utf-8") as f:
        f.write("".join(morceaux))
    return chemin
