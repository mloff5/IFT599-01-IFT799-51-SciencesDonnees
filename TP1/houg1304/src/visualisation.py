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
# Fichier : src/visualisation.py
# Role    : production des figures de la Methode 2 (nuages de points avec
#           distributions marginales).
#
# Sources : matplotlib (trace des figures). La mise en page de type "jointplot"
#           (nuage de points central et histogrammes marginaux) est reproduite
#           avec matplotlib.gridspec par l'equipe, en s'inspirant de la figure
#           de reference citee dans l'enonce :
#           https://seaborn.pydata.org/generated/seaborn.jointplot.html
# =============================================================================
"""Figures : nuages de points colories par classe et histogrammes marginaux."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # rendu hors ecran : aucune fenetre n'est ouverte
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config  # noqa: E402


def _couleur(classe, index=0):
    palette_secours = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd"]
    return config.COULEURS_CLASSES.get(classe, palette_secours[index % len(palette_secours)])


def _marqueur(classe, index=0):
    marqueurs_secours = ["o", "s", "^", "D", "v"]
    return config.MARQUEURS_CLASSES.get(classe, marqueurs_secours[index % len(marqueurs_secours)])


def nuage_de_points(coordonnees, y, classes, titre, etiquette_x, etiquette_y,
                    chemin, sous_titre=None, histogrammes=None):
    """Trace un nuage de points a deux variables, colorie par classe.

    Si `histogrammes` est vrai, les distributions marginales de chaque classe
    sont affichees le long des deux axes (point 2.b optionnel de l'enonce :
    illustration de la separation des classes par un affichage 1D).

    Parametres
    ----------
    coordonnees : ndarray (m, 2)
    y : ndarray (m,) etiquettes de classe
    classes : liste ordonnee des classes a afficher
    chemin : fichier de sortie (png)
    """
    if histogrammes is None:
        histogrammes = config.AFFICHER_HISTOGRAMMES

    coordonnees = np.asarray(coordonnees, dtype=np.float64)
    y = np.asarray(y)
    chemin = Path(chemin)
    chemin.parent.mkdir(parents=True, exist_ok=True)

    if histogrammes:
        figure = plt.figure(figsize=(8.0, 7.4))
        grille = figure.add_gridspec(
            2, 2, width_ratios=(5, 1.4), height_ratios=(1.4, 5),
            left=0.10, right=0.97, bottom=0.09, top=0.88,
            wspace=0.05, hspace=0.05,
        )
        axe = figure.add_subplot(grille[1, 0])
        axe_haut = figure.add_subplot(grille[0, 0], sharex=axe)
        axe_droite = figure.add_subplot(grille[1, 1], sharey=axe)
        axe_haut.tick_params(axis="x", labelbottom=False)
        axe_droite.tick_params(axis="y", labelleft=False)
        for a in (axe_haut, axe_droite):
            a.spines["top"].set_visible(False)
            a.spines["right"].set_visible(False)
    else:
        figure, axe = plt.subplots(figsize=(7.6, 6.4))
        axe_haut = axe_droite = None

    for index, classe in enumerate(classes):
        masque = (y == classe)
        if not np.any(masque):
            continue
        couleur = _couleur(classe, index)
        axe.scatter(
            coordonnees[masque, 0], coordonnees[masque, 1],
            s=26, alpha=0.80, c=couleur, marker=_marqueur(classe, index),
            edgecolors="white", linewidths=0.4,
            label="{} (n={})".format(classe, int(masque.sum())),
        )
        if histogrammes:
            n_bacs = max(10, min(30, int(np.sqrt(masque.sum()) * 2)))
            axe_haut.hist(coordonnees[masque, 0], bins=n_bacs, color=couleur,
                          alpha=0.55, density=True, edgecolor="none")
            axe_droite.hist(coordonnees[masque, 1], bins=n_bacs, color=couleur,
                            alpha=0.55, density=True, edgecolor="none",
                            orientation="horizontal")

    axe.set_xlabel(etiquette_x)
    axe.set_ylabel(etiquette_y)
    axe.grid(True, linestyle=":", linewidth=0.6, alpha=0.6)
    axe.set_axisbelow(True)
    axe.legend(loc="best", frameon=True, framealpha=0.9, fontsize=9)

    if sous_titre:
        figure.suptitle(titre, fontsize=12.5, fontweight="bold", y=0.975)
        cible = axe_haut if histogrammes else axe
        cible.set_title(sous_titre, fontsize=9.5, color="#333333", pad=6)
    else:
        figure.suptitle(titre, fontsize=12.5, fontweight="bold", y=0.975)

    figure.savefig(chemin, dpi=config.DPI_FIGURES, bbox_inches="tight")
    plt.close(figure)
    return chemin


def histogrammes_1d(valeurs_par_classe, titre, etiquette_x, chemin):
    """Affichage 1D : distribution de chaque classe pour une seule variable."""
    chemin = Path(chemin)
    chemin.parent.mkdir(parents=True, exist_ok=True)
    figure, axe = plt.subplots(figsize=(7.4, 4.2))
    for index, (classe, valeurs) in enumerate(valeurs_par_classe.items()):
        n_bacs = max(10, min(30, int(np.sqrt(len(valeurs)) * 2)))
        axe.hist(valeurs, bins=n_bacs, alpha=0.55, density=True,
                 color=_couleur(classe, index), label=classe, edgecolor="none")
    axe.set_xlabel(etiquette_x)
    axe.set_ylabel("Densite")
    axe.set_title(titre, fontsize=11.5, fontweight="bold")
    axe.grid(True, linestyle=":", linewidth=0.6, alpha=0.6)
    axe.set_axisbelow(True)
    axe.legend(frameon=True, fontsize=9)
    figure.tight_layout()
    figure.savefig(chemin, dpi=config.DPI_FIGURES)
    plt.close(figure)
    return chemin


def courbe_variance_acp(variance_expliquee, variance_cumulee, chemin):
    """Eboulis des valeurs propres de l'ACP (figure d'appui pour le rapport)."""
    chemin = Path(chemin)
    chemin.parent.mkdir(parents=True, exist_ok=True)
    composantes = np.arange(1, len(variance_expliquee) + 1)

    figure, axe = plt.subplots(figsize=(7.4, 4.4))
    axe.bar(composantes, variance_expliquee * 100.0, color="#4c78a8",
            alpha=0.85, label="Variance expliquee par composante")
    axe.set_xlabel("Composante principale")
    axe.set_ylabel("Variance expliquee (%)")
    axe.set_xticks(composantes)
    axe.grid(True, axis="y", linestyle=":", linewidth=0.6, alpha=0.6)
    axe.set_axisbelow(True)

    axe_cumul = axe.twinx()
    axe_cumul.plot(composantes, variance_cumulee * 100.0, color="#d62728",
                   marker="o", markersize=4, linewidth=1.6,
                   label="Variance cumulee")
    axe_cumul.set_ylabel("Variance cumulee (%)")
    axe_cumul.set_ylim(0, 100)

    lignes_1, etiquettes_1 = axe.get_legend_handles_labels()
    lignes_2, etiquettes_2 = axe_cumul.get_legend_handles_labels()
    axe.legend(lignes_1 + lignes_2, etiquettes_1 + etiquettes_2,
               loc="center right", fontsize=9)
    axe.set_title("ACP - variance expliquee par les {} premieres composantes"
                  .format(len(variance_expliquee)),
                  fontsize=11.5, fontweight="bold")
    figure.tight_layout()
    figure.savefig(chemin, dpi=config.DPI_FIGURES)
    plt.close(figure)
    return chemin


def carte_de_chaleur_overlap(matrice, classes, titre, chemin, seuil=1.0):
    """Carte de chaleur des valeurs de Overlap entre toutes les paires."""
    chemin = Path(chemin)
    chemin.parent.mkdir(parents=True, exist_ok=True)
    matrice = np.asarray(matrice, dtype=np.float64)

    figure, axe = plt.subplots(figsize=(6.6, 5.6))
    donnees = np.ma.masked_invalid(matrice)
    image = axe.imshow(donnees, cmap="RdYlGn_r", vmin=0.0,
                       vmax=max(2.0, float(np.nanmax(matrice[np.isfinite(matrice)]))
                                if np.any(np.isfinite(matrice)) else 2.0))
    axe.set_xticks(range(len(classes)))
    axe.set_yticks(range(len(classes)))
    axe.set_xticklabels(classes)
    axe.set_yticklabels(classes)

    for i in range(len(classes)):
        for j in range(len(classes)):
            if i == j:
                axe.text(j, i, "-", ha="center", va="center", fontsize=9, color="#555555")
                continue
            valeur = matrice[i, j]
            if not np.isfinite(valeur):
                texte = "inf"
            else:
                texte = "{:.2f}".format(valeur)
            axe.text(j, i, texte, ha="center", va="center", fontsize=9,
                     fontweight="bold" if np.isfinite(valeur) and valeur < seuil else "normal",
                     color="black")

    axe.set_title(titre, fontsize=11.5, fontweight="bold")
    figure.colorbar(image, ax=axe, shrink=0.85, label="Overlap")
    figure.tight_layout()
    figure.savefig(chemin, dpi=config.DPI_FIGURES)
    plt.close(figure)
    return chemin
