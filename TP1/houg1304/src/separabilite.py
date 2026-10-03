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
# Fichier : src/separabilite.py
# Role    : mesures de separabilite de classes (cohesion, separation, Overlap)
#           telles que definies a la section 1.2 de l'enonce.
#
# Sources : numpy uniquement. Les formules sont celles de l'enonce du TP1.
# =============================================================================
"""Cohesion, separation et indicateur de superposition Overlap().

Definitions (enonce du TP1, section 1.2)
----------------------------------------
Soit une classe C = {x_1, ..., x_n} de centre x_C (moyenne des points de C).

* Distance intra-classe (cohesion) :
      dist_intra(C) = max { dist(x_i, x_C) | x_i in C }

* Distance inter-classe (separation) entre C1 = {x_1..x_n1} et C2 = {y_1..y_n2} :
      dist(C1, C2)  = min { dist(x_i, y_C2) | x_i in C1 }
      dist(C2, C1)  = min { dist(y_j, x_C1) | y_j in C2 }
      dist_inter(C1, C2) = min( dist(C1, C2), dist(C2, C1) )

* Indicateur de superposition :
      Overlap(C1, C2) = (dist_intra(C1) + dist_intra(C2)) / (2 * dist_inter(C1, C2))

  Overlap(C1, C2) < 1  =>  les classes C1 et C2 sont bien separees.
  (condition suffisante, non necessaire)

Toutes les fonctions de ce module travaillent dans un espace ou la distance a
utiliser est la distance euclidienne : pour la distance de Mahalanobis, les
donnees sont prealablement blanchies par `distances.preparer_espace`.
"""

from __future__ import annotations

import numpy as np

from .distances import distances_au_point


# -----------------------------------------------------------------------------
# Mesures elementaires
# -----------------------------------------------------------------------------
def centre_classe(Xc):
    """Centre x_C d'une classe : moyenne des points de la classe."""
    return np.asarray(Xc, dtype=np.float64).mean(axis=0)


def distance_intra(Xc, centre=None):
    """dist_intra(C) = max distance d'un point de C au centre de C."""
    Xc = np.asarray(Xc, dtype=np.float64)
    if Xc.shape[0] == 0:
        return np.nan
    if centre is None:
        centre = centre_classe(Xc)
    return float(distances_au_point(Xc, centre).max())


def distance_inter(X1, X2, centre1=None, centre2=None):
    """dist_inter(C1, C2) = min( dist(C1,C2), dist(C2,C1) ).

    dist(C1, C2) est la distance minimale entre un objet quelconque de C1 et le
    centre de C2 ; dist(C2, C1) est definie symetriquement.
    """
    X1 = np.asarray(X1, dtype=np.float64)
    X2 = np.asarray(X2, dtype=np.float64)
    if centre1 is None:
        centre1 = centre_classe(X1)
    if centre2 is None:
        centre2 = centre_classe(X2)
    d_12 = float(distances_au_point(X1, centre2).min())
    d_21 = float(distances_au_point(X2, centre1).min())
    return min(d_12, d_21)


def overlap(X1, X2):
    """Indicateur de superposition Overlap(C1, C2)."""
    centre1 = centre_classe(X1)
    centre2 = centre_classe(X2)
    intra1 = distance_intra(X1, centre1)
    intra2 = distance_intra(X2, centre2)
    inter = distance_inter(X1, X2, centre1, centre2)
    if inter <= 0.0:
        return np.inf
    return (intra1 + intra2) / (2.0 * inter)


def details_paire(X1, X2):
    """Retourne toutes les quantites intermediaires pour une paire de classes."""
    centre1 = centre_classe(X1)
    centre2 = centre_classe(X2)
    intra1 = distance_intra(X1, centre1)
    intra2 = distance_intra(X2, centre2)
    d_12 = float(distances_au_point(X1, centre2).min())
    d_21 = float(distances_au_point(X2, centre1).min())
    inter = min(d_12, d_21)
    valeur = np.inf if inter <= 0.0 else (intra1 + intra2) / (2.0 * inter)
    return {
        "dist_intra_C1": intra1,
        "dist_intra_C2": intra2,
        "dist_C1_vers_C2": d_12,
        "dist_C2_vers_C1": d_21,
        "dist_inter": inter,
        "overlap": valeur,
        "bien_separees": bool(valeur < 1.0),
    }


# -----------------------------------------------------------------------------
# Calculs sur l'ensemble des classes
# -----------------------------------------------------------------------------
def toutes_les_cohesions(Z, y, classes):
    """dist_intra pour chaque classe."""
    resultat = {}
    for c in classes:
        Xc = Z[y == c]
        resultat[c] = distance_intra(Xc)
    return resultat


def toutes_les_paires(Z, y, classes):
    """Calcule les details de separabilite pour toutes les paires mutuelles.

    Retour
    ------
    liste de dictionnaires, un par paire de classes, triee par classe.
    """
    classes = sorted(classes)
    sous_ensembles = {c: Z[y == c] for c in classes}
    resultats = []
    for i in range(len(classes)):
        for j in range(i + 1, len(classes)):
            c1, c2 = classes[i], classes[j]
            details = details_paire(sous_ensembles[c1], sous_ensembles[c2])
            details["classe_1"] = c1
            details["classe_2"] = c2
            resultats.append(details)
    return resultats


# -----------------------------------------------------------------------------
# Versions vectorisees pour la selection de variables (Methode 2)
# -----------------------------------------------------------------------------
def overlap_1d_toutes_variables(X, indices_c1, indices_c2):
    """Overlap univarie calcule simultanement pour toutes les variables.

    Cette version vectorisee evalue, en une seule passe, la mesure Overlap()
    appliquee a chacune des p variables prise isolement (espace a 1 dimension).
    Elle sert d'etape de filtrage dans la strategie de recherche de la Methode 2.

    Parametres
    ----------
    X : ndarray (n, p)
    indices_c1, indices_c2 : indices des lignes des deux classes

    Retour
    ------
    ndarray (p,) des valeurs de Overlap, variable par variable.
    """
    X1 = X[indices_c1]
    X2 = X[indices_c2]
    centre1 = X1.mean(axis=0)
    centre2 = X2.mean(axis=0)

    intra1 = np.abs(X1 - centre1[None, :]).max(axis=0)
    intra2 = np.abs(X2 - centre2[None, :]).max(axis=0)

    d_12 = np.abs(X1 - centre2[None, :]).min(axis=0)
    d_21 = np.abs(X2 - centre1[None, :]).min(axis=0)
    inter = np.minimum(d_12, d_21)

    with np.errstate(divide="ignore", invalid="ignore"):
        valeurs = (intra1 + intra2) / (2.0 * inter)
    valeurs[~np.isfinite(valeurs)] = np.inf
    return valeurs


def overlap_2d(X, j, k, indices_c1, indices_c2):
    """Overlap calcule sur le couple de variables (j, k) uniquement."""
    colonnes = [j, k]
    return overlap(X[np.ix_(indices_c1, colonnes)],
                   X[np.ix_(indices_c2, colonnes)])


def overlap_2d_colonne_fixee(X, k, indices_c1, indices_c2):
    """Overlap de toutes les paires (j, k), la variable k etant fixee.

    Cette version vectorisee evalue en une seule passe les p paires de
    variables contenant la variable k. Elle repose sur le fait que la racine
    carree est croissante : pour un maximum ou un minimum, il suffit de
    travailler sur les carres des distances.

        dist((x_ij, x_ik), (c_j, c_k))^2 = (x_ij - c_j)^2 + (x_ik - c_k)^2

    Le second terme ne depend pas de j : il se diffuse par broadcasting sur
    toutes les colonnes. Le cout d'une passe est donc celui de quelques
    operations sur des matrices (n, p), ce qui permet a la recherche locale de
    balayer l'integralite des variables a chaque iteration au lieu de se
    limiter a la liste de candidats.

    Parametres
    ----------
    X : ndarray (n, p)
    k : index de la variable fixee
    indices_c1, indices_c2 : indices des lignes des deux classes

    Retour
    ------
    ndarray (p,) ; la case d'index k vaut l'infini (paire degeneree j = k).
    """
    X1 = X[indices_c1]
    X2 = X[indices_c2]
    centre1 = X1.mean(axis=0)
    centre2 = X2.mean(axis=0)

    # Cohesion : max sur les objets de la classe, pour chaque variable j.
    ecarts1 = X1 - centre1[None, :]
    ecarts2 = X2 - centre2[None, :]
    intra1 = np.sqrt((ecarts1 ** 2 + (ecarts1[:, k] ** 2)[:, None]).max(axis=0))
    intra2 = np.sqrt((ecarts2 ** 2 + (ecarts2[:, k] ** 2)[:, None]).max(axis=0))

    # Separation : min sur les objets d'une classe vers le centre de l'autre.
    croises1 = X1 - centre2[None, :]
    croises2 = X2 - centre1[None, :]
    d_12 = np.sqrt((croises1 ** 2 + (croises1[:, k] ** 2)[:, None]).min(axis=0))
    d_21 = np.sqrt((croises2 ** 2 + (croises2[:, k] ** 2)[:, None]).min(axis=0))
    inter = np.minimum(d_12, d_21)

    with np.errstate(divide="ignore", invalid="ignore"):
        valeurs = (intra1 + intra2) / (2.0 * inter)
    valeurs[~np.isfinite(valeurs)] = np.inf
    valeurs[k] = np.inf  # une paire doit etre formee de deux variables distinctes
    return valeurs
