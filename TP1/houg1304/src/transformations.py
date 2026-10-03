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
# Fichier : src/transformations.py
# Role    : transformations des donnees (min-max, Z-score), ACP et t-SNE.
#
# Sources :
#   * Min-max et Z-score : diapositives 34 et 35 du Theme 1 du cours.
#   * ACP : implantation de l'equipe avec numpy, suivant exactement la procedure
#     des diapositives 38 et 39 du Theme 1 (centrage-reduction, matrice de
#     correlation, valeurs/vecteurs propres, projection sur les axes
#     principaux). Aucune bibliotheque externe d'ACP n'est utilisee.
#   * t-SNE : scikit-learn, classe sklearn.manifold.TSNE.
#     Reference : van der Maaten, L. & Hinton, G. (2008), "Visualizing Data
#     using t-SNE", Journal of Machine Learning Research 9:2579-2605.
#     Documentation : https://scikit-learn.org/stable/modules/generated/
#                     sklearn.manifold.TSNE.html
# =============================================================================
"""Transformations de contenu et reductions de dimension."""

from __future__ import annotations

import numpy as np


# -----------------------------------------------------------------------------
# Transformations de contenu (Theme 1, diapositives 34-35)
# -----------------------------------------------------------------------------
def brut(X):
    """Aucune transformation : les donnees sont utilisees telles quelles."""
    return np.asarray(X, dtype=np.float64)


def zscore(X):
    """Standardisation Z-score : x' = (x - mu) / sigma, colonne par colonne."""
    X = np.asarray(X, dtype=np.float64)
    moyenne = X.mean(axis=0, keepdims=True)
    ecart = X.std(axis=0, keepdims=True)
    ecart = np.where(ecart > 0.0, ecart, 1.0)  # colonne constante -> inchangee
    return (X - moyenne) / ecart


def minmax(X, nouveau_min=0.0, nouveau_max=1.0):
    """Normalisation min-max dans l'intervalle [nouveau_min, nouveau_max]."""
    X = np.asarray(X, dtype=np.float64)
    ancien_min = X.min(axis=0, keepdims=True)
    ancien_max = X.max(axis=0, keepdims=True)
    etendue = ancien_max - ancien_min
    etendue = np.where(etendue > 0.0, etendue, 1.0)
    return (X - ancien_min) / etendue * (nouveau_max - nouveau_min) + nouveau_min


TRANSFORMATIONS = {
    "brut": brut,
    "zscore": zscore,
    "minmax": minmax,
}

ETIQUETTES_TRANSFORMATIONS = {
    "brut": "Aucune transformation",
    "zscore": "Standardisation Z-score",
    "minmax": "Normalisation min-max [0,1]",
}


def appliquer_transformation(X, nom):
    if nom not in TRANSFORMATIONS:
        raise ValueError("Transformation inconnue : {}".format(nom))
    return TRANSFORMATIONS[nom](X)


# -----------------------------------------------------------------------------
# Analyse en composantes principales (Theme 1, diapositives 38-39)
# -----------------------------------------------------------------------------
class ACP:
    """Analyse en composantes principales implantee avec numpy.

    Procedure suivie (diapositives 38-39 du Theme 1) :
      1. Centrage-reduction de la matrice A -> A' (methode Z-score par colonne).
      2. Calcul de la matrice de correlation V = (A')^T A' / (n-1).
      3. Calcul des valeurs propres lambda_k et vecteurs propres u_k de V.
      4. Conservation des vecteurs propres associes aux plus grandes valeurs
         propres (axes principaux).
      5. Projection des objets sur chaque axe principal : pi_k(M) = M . u_k.

    Le calcul est effectue par decomposition en valeurs singulieres (SVD) de A',
    numeriquement equivalente a la diagonalisation de V mais plus stable et plus
    rapide lorsque p > n (ici 2000 variables pour 801 patients).
    """

    def __init__(self, n_composantes=15, centrer_reduire=True):
        self.n_composantes = n_composantes
        self.centrer_reduire = centrer_reduire
        self.moyennes_ = None
        self.ecarts_ = None
        self.axes_ = None                 # (n_composantes, p) : vecteurs propres u_k
        self.valeurs_propres_ = None      # (n_composantes,)   : lambda_k
        self.variance_expliquee_ = None   # proportion par composante
        self.variance_cumulee_ = None

    def ajuster(self, X):
        X = np.asarray(X, dtype=np.float64)
        n = X.shape[0]

        # Etape 1 : centrage (et reduction si demande).
        self.moyennes_ = X.mean(axis=0, keepdims=True)
        A = X - self.moyennes_
        if self.centrer_reduire:
            ecarts = X.std(axis=0, keepdims=True)
            self.ecarts_ = np.where(ecarts > 0.0, ecarts, 1.0)
            A = A / self.ecarts_
        else:
            self.ecarts_ = np.ones_like(self.moyennes_)

        # Etapes 2 a 4 : valeurs et vecteurs propres de V = A^T A / (n-1),
        # obtenus par SVD de A = U S Vt  =>  lambda_k = s_k^2 / (n-1), u_k = Vt_k.
        _, valeurs_singulieres, Vt = np.linalg.svd(A, full_matrices=False)
        valeurs_propres = (valeurs_singulieres ** 2) / max(n - 1, 1)

        k = min(self.n_composantes, Vt.shape[0])
        self.axes_ = Vt[:k]
        self.valeurs_propres_ = valeurs_propres[:k]

        variance_totale = valeurs_propres.sum()
        if variance_totale > 0:
            self.variance_expliquee_ = self.valeurs_propres_ / variance_totale
        else:
            self.variance_expliquee_ = np.zeros(k)
        self.variance_cumulee_ = np.cumsum(self.variance_expliquee_)
        return self

    def transformer(self, X):
        """Etape 5 : projection des objets sur les axes principaux."""
        X = np.asarray(X, dtype=np.float64)
        A = (X - self.moyennes_) / self.ecarts_
        return A @ self.axes_.T

    def ajuster_transformer(self, X):
        return self.ajuster(X).transformer(X)

    def noms_composantes(self):
        return ["CP{}".format(i + 1) for i in range(self.axes_.shape[0])]


# -----------------------------------------------------------------------------
# t-SNE (scikit-learn)
# -----------------------------------------------------------------------------
def tsne(X, n_composantes=2, perplexite=30.0, graine=42,
         n_iter=1000, learning_rate="auto", verbeux=False):
    """Projection t-SNE en `n_composantes` dimensions.

    Source : scikit-learn, sklearn.manifold.TSNE (van der Maaten & Hinton, 2008).
    """
    from sklearn.manifold import TSNE

    X = np.asarray(X, dtype=np.float64)
    n = X.shape[0]
    # La perplexite doit rester strictement inferieure au nombre d'objets.
    perplexite_effective = float(min(perplexite, max(5.0, (n - 1) / 3.0)))

    parametres = dict(
        n_components=n_composantes,
        perplexity=perplexite_effective,
        learning_rate=learning_rate,
        init="pca",
        random_state=graine,
        verbose=1 if verbeux else 0,
    )
    # Le nom du parametre du nombre d'iterations a change selon les versions de
    # scikit-learn (n_iter -> max_iter). On s'adapte automatiquement.
    try:
        modele = TSNE(max_iter=n_iter, **parametres)
    except TypeError:
        modele = TSNE(n_iter=n_iter, **parametres)
    return modele.fit_transform(X), perplexite_effective
