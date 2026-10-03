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
# Fichier : src/distances.py
# Role    : mesures de distance (Euclidienne et Mahalanobis).
#
# Sources : numpy uniquement. L'implantation de Mahalanobis par blanchiment
#           (whitening) est ecrite par l'equipe a partir de la definition
#           d_M(x,y) = sqrt((x-y)^T Sigma^-1 (x-y)).
# =============================================================================
"""Mesures de distance utilisees par les deux methodes du TP.

Idee centrale
-------------
La distance de Mahalanobis est une distance euclidienne apres une
transformation lineaire des donnees. En effet, si Sigma^-1 = W^T W alors

    d_M(x, y) = sqrt((x-y)^T W^T W (x-y)) = || W x - W y ||_2

Il suffit donc de calculer une seule fois la matrice de blanchiment W, de
transformer la matrice de donnees (Z = X W^T), puis de n'utiliser que la
distance euclidienne dans tout le reste du programme. Cela evite de recalculer
une forme quadratique pour chaque couple de points, ce qui serait prohibitif
avec 801 patients et 2000 variables.

Remarque importante sur la covariance
-------------------------------------
Avec n = 801 patients et p = 2000 variables, la matrice de covariance empirique
est de rang au plus n-1 = 800 : elle n'est pas inversible. Deux strategies sont
proposees (voir config.METHODE_INVERSION_COVARIANCE) :

* "ridge" : on regularise Sigma_reg = Sigma + lambda * (trace(Sigma)/p) * I.
            La matrice devient definie positive et le conditionnement est
            controle par lambda.
* "pinv"  : pseudo-inverse de Moore-Penrose. Le blanchiment est alors exact sur
            le sous-espace de rang plein, mais toutes les directions de ce
            sous-espace sont ramenees a une variance unitaire, ce qui peut
            uniformiser artificiellement les distances.
"""

from __future__ import annotations

import numpy as np


# -----------------------------------------------------------------------------
# Distance euclidienne
# -----------------------------------------------------------------------------
def distances_au_point(A, v):
    """Distances euclidiennes entre chaque ligne de A et le vecteur v.

    Parametres
    ----------
    A : ndarray (m, p)
    v : ndarray (p,)

    Retour
    ------
    ndarray (m,) des distances euclidiennes.
    """
    difference = A - np.asarray(v, dtype=np.float64).reshape(1, -1)
    return np.sqrt(np.einsum("ij,ij->i", difference, difference))


def matrice_distances(A, B=None):
    """Matrice des distances euclidiennes entre les lignes de A et celles de B."""
    A = np.asarray(A, dtype=np.float64)
    B = A if B is None else np.asarray(B, dtype=np.float64)
    carres_A = np.einsum("ij,ij->i", A, A)[:, None]
    carres_B = np.einsum("ij,ij->i", B, B)[None, :]
    produit = A @ B.T
    carres = np.maximum(carres_A + carres_B - 2.0 * produit, 0.0)
    return np.sqrt(carres)


# -----------------------------------------------------------------------------
# Distance de Mahalanobis (par blanchiment)
# -----------------------------------------------------------------------------
def matrice_blanchiment(X, methode="ridge", lambda_ridge=1e-2):
    """Calcule la matrice de blanchiment W telle que d_M(x,y) = ||Wx - Wy||.

    La covariance utilisee est la covariance globale du jeu de donnees (toutes
    classes confondues). C'est le choix usuel lorsque l'on veut comparer des
    classes entre elles avec une metrique commune : utiliser une covariance
    propre a chaque classe rendrait les distances non comparables d'une classe
    a l'autre (la mesure ne serait plus une distance unique).

    Parametres
    ----------
    X : ndarray (n, p)
    methode : "ridge" ou "pinv"
    lambda_ridge : intensite de la regularisation (methode "ridge")

    Retour
    ------
    W : ndarray (k, p) telle que Z = X @ W.T donne un espace ou la distance
        euclidienne vaut la distance de Mahalanobis dans l'espace d'origine.
    infos : dictionnaire d'informations de diagnostic.
    """
    X = np.asarray(X, dtype=np.float64)
    n, p = X.shape
    centre = X - X.mean(axis=0, keepdims=True)
    # Covariance empirique non biaisee.
    sigma = (centre.T @ centre) / max(n - 1, 1)

    # Decomposition spectrale (sigma est symetrique semi-definie positive).
    valeurs, vecteurs = np.linalg.eigh(sigma)
    valeurs = np.maximum(valeurs, 0.0)
    rang_numerique = int(np.sum(valeurs > valeurs.max() * 1e-12)) if valeurs.max() > 0 else 0

    if methode == "ridge":
        delta = lambda_ridge * (valeurs.sum() / p) if valeurs.sum() > 0 else lambda_ridge
        valeurs_reg = valeurs + delta
        inverses = 1.0 / np.sqrt(valeurs_reg)
        W = (vecteurs * inverses[None, :]).T
        infos = {
            "methode": "ridge",
            "lambda": lambda_ridge,
            "delta_absolu": float(delta),
            "rang_numerique": rang_numerique,
            "dimension_blanchie": int(W.shape[0]),
            "conditionnement": float(valeurs_reg.max() / valeurs_reg.min()),
        }
    elif methode == "pinv":
        seuil = valeurs.max() * 1e-12 if valeurs.max() > 0 else 0.0
        garder = valeurs > seuil
        inverses = 1.0 / np.sqrt(valeurs[garder])
        W = (vecteurs[:, garder] * inverses[None, :]).T
        infos = {
            "methode": "pinv",
            "lambda": 0.0,
            "delta_absolu": 0.0,
            "rang_numerique": rang_numerique,
            "dimension_blanchie": int(W.shape[0]),
            "conditionnement": float(valeurs[garder].max() / valeurs[garder].min()),
        }
    else:
        raise ValueError("Methode d'inversion inconnue : {}".format(methode))

    return W, infos


def appliquer_blanchiment(X, W):
    """Projette X dans l'espace blanchi : Z = X @ W.T."""
    return np.asarray(X, dtype=np.float64) @ np.asarray(W, dtype=np.float64).T


def preparer_espace(X, metrique, methode_inversion="ridge", lambda_ridge=1e-2):
    """Retourne la representation de X dans laquelle la distance euclidienne
    coincide avec la metrique demandee.

    Parametres
    ----------
    metrique : "euclidienne" ou "mahalanobis"

    Retour
    ------
    (Z, infos) ou Z est la matrice a utiliser avec `distances_au_point`.
    """
    if metrique == "euclidienne":
        return np.asarray(X, dtype=np.float64), {"methode": "aucune"}
    if metrique == "mahalanobis":
        W, infos = matrice_blanchiment(X, methode_inversion, lambda_ridge)
        return appliquer_blanchiment(X, W), infos
    raise ValueError("Metrique inconnue : {}".format(metrique))
