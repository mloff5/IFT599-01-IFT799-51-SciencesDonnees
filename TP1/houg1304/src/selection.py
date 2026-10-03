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
# Fichier : src/selection.py
# Role    : selection automatique de paires de variables (Methode 2).
#
# Sources : numpy uniquement. L'algorithme suit le schema general de selection
#           de dimensions presente a la diapositive 37 du Theme 1 du cours ;
#           la strategie de recherche et les mesures d'evaluation ont ete
#           definies par l'equipe.
# =============================================================================
"""Selection automatique de la meilleure paire de variables.

Schema implante (diapositive 37 du Theme 1)
-------------------------------------------

    Attributs --> Strategie de recherche --> Sous-ensemble d'attributs
                        ^                            |
                        |                            v
                  Critere d'arret <------------- Evaluation
                        | oui
                        v
                 Attributs selectionnes --> Processus de validation

Les quatre briques sont materialisees ainsi :

* Attributs             : les p variables de l'ensemble de depart.
* Strategie de recherche: en trois temps.
                          (A) filtrage univarie : la mesure d'evaluation est
                              appliquee a chaque variable prise seule, et les
                              N_CANDIDATES meilleures sont retenues ;
                          (B) pour chaque variable candidate k, une passe
                              vectorisee evalue toutes les paires (j, k) avec
                              j parcourant les p variables ;
                          (C) recherche locale par echange de coordonnees : on
                              fixe une variable et on cherche la meilleure
                              partenaire parmi les p variables, en alternant,
                              jusqu'au critere d'arret.
                          Avec p = 2000 et N_CANDIDATES = 60, environ 6 % des
                          p(p-1)/2 paires possibles sont examinees, pour un
                          cout de quelques dizaines d'operations matricielles.
                          On reste donc tres loin d'une recherche exhaustive,
                          conformement a la note 2 de l'enonce, tout en evitant
                          les mauvais optima locaux d'une recherche gloutonne
                          confinee a la liste de candidats.
* Evaluation            : une mesure construite a partir de Overlap() (deux
                          mesures distinctes, voir ci-dessous).
* Critere d'arret       : plus aucune amelioration relative superieure a
                          EPSILON, ou nombre maximal d'iterations atteint.
* Validation            : verification de la condition Overlap < 1 et
                          production de la figure correspondante.

Deux mesures d'evaluation distinctes
------------------------------------
L'enonce precise que les mesures d'evaluation ne doivent pas etre les memes
selon l'objectif poursuivi :

1. Meilleure paire pour UNE paire de classes donnee (C_a, C_b) :
       f(j,k) = Overlap_{a,b}(j,k)                       -> a minimiser.

2. Meilleure paire pour TOUTES les paires de classes simultanement. Deux
   mesures agregees, toutes deux construites a partir de Overlap(), sont
   implantees et comparees :

       g_max(j,k)     = max over (a,b) of Overlap_{a,b}(j,k)   -> a minimiser
       g_moyenne(j,k) = mean over (a,b) of Overlap_{a,b}(j,k)  -> a minimiser

   * g_max est un critere minimax (pire cas) : une paire de variables n'est
     jugee bonne que si elle separe convenablement meme les deux classes les
     plus difficiles a distinguer. C'est le critere le plus exigeant, mais il
     est tres sensible a une seule paire de classes pathologique.
   * g_moyenne privilegie la qualite d'ensemble de l'affichage : elle accepte
     qu'une paire de classes reste confondue si toutes les autres sont bien
     separees. C'est generalement ce critere qui produit la figure la plus
     lisible pour les cinq classes affichees simultanement.

   Les deux criteres sont executes et leurs resultats compares dans le rapport.
"""

from __future__ import annotations

import numpy as np

from .separabilite import (
    overlap_1d_toutes_variables, overlap_2d, overlap_2d_colonne_fixee,
)


# -----------------------------------------------------------------------------
# Mesures d'evaluation
# -----------------------------------------------------------------------------
class _EvaluateurBase:
    """Comportement commun aux mesures d'evaluation.

    Trois modes d'evaluation sont exposes :

    * `evaluer_1d(X)`              -> score de chaque variable prise seule ;
    * `evaluer_2d(X, j, k)`        -> score d'une paire de variables ;
    * `evaluer_colonne_fixee(X, k)`-> score de toutes les paires (j, k), la
      variable k etant fixee, en une seule passe vectorisee.

    Le troisieme mode est celui qui rend la recherche efficace : il permet de
    balayer l'integralite des p variables au cout de quelques operations
    matricielles, au lieu de p appels successifs.
    """

    def candidats(self, X, k):
        scores = self.evaluer_1d(X)
        k = int(min(k, X.shape[1]))
        return [int(i) for i in np.argsort(scores, kind="stable")[:k]]


class EvaluateurPaireDeClasses(_EvaluateurBase):
    """Mesure d'evaluation f(j,k) = Overlap(C_a, C_b) sur les variables (j,k)."""

    def __init__(self, y, classe_1, classe_2):
        self.classe_1 = classe_1
        self.classe_2 = classe_2
        self.indices_1 = np.flatnonzero(y == classe_1)
        self.indices_2 = np.flatnonzero(y == classe_2)
        self.nom = "Overlap({}, {})".format(classe_1, classe_2)
        self.cible = "{}-{}".format(classe_1, classe_2)

    def evaluer_1d(self, X):
        """Vecteur des scores univaries (une valeur par variable)."""
        return overlap_1d_toutes_variables(X, self.indices_1, self.indices_2)

    def evaluer_2d(self, X, j, k):
        return overlap_2d(X, j, k, self.indices_1, self.indices_2)

    def evaluer_colonne_fixee(self, X, k):
        return overlap_2d_colonne_fixee(X, k, self.indices_1, self.indices_2)

    def detail(self, X, j, k):
        return {"overlap": self.evaluer_2d(X, j, k)}


class EvaluateurGlobal(_EvaluateurBase):
    """Mesure d'evaluation agregee sur l'ensemble des paires de classes.

    Parametres
    ----------
    agregation : "max" (critere minimax, pire cas) ou "moyenne".
    """

    AGREGATIONS = {
        "max": ("max des Overlap sur les {} paires de classes", np.maximum, max),
        "moyenne": ("moyenne des Overlap sur les {} paires de classes", None, None),
    }

    def __init__(self, y, classes, agregation="max"):
        if agregation not in self.AGREGATIONS:
            raise ValueError("Agregation inconnue : {}".format(agregation))
        self.agregation = agregation
        self.classes = sorted(classes)
        self.paires = [(self.classes[i], self.classes[j])
                       for i in range(len(self.classes))
                       for j in range(i + 1, len(self.classes))]
        self.indices = {c: np.flatnonzero(y == c) for c in self.classes}
        self.nom = self.AGREGATIONS[agregation][0].format(len(self.paires))
        self.cible = "toutes"

    def _agreger(self, valeurs):
        valeurs = np.asarray(valeurs, dtype=np.float64)
        if self.agregation == "max":
            return valeurs.max(axis=0)
        return valeurs.mean(axis=0)

    def _scores_1d_par_paire(self, X):
        return np.vstack([
            overlap_1d_toutes_variables(X, self.indices[c1], self.indices[c2])
            for c1, c2 in self.paires
        ])

    def evaluer_1d(self, X):
        return self._agreger(self._scores_1d_par_paire(X))

    def candidats(self, X, k):
        """Sous-ensemble d'attributs examine par la recherche locale.

        Un filtrage fonde uniquement sur le score agrege est trompeur : le
        classement est alors entierement dicte par la paire de classes la plus
        difficile a separer, et les variables utiles aux autres paires
        disparaissent du sous-ensemble. La liste de candidats est donc
        diversifiee :

          * un tiers du budget est attribue aux meilleures variables au sens du
            score agrege ;
          * les deux tiers restants sont repartis equitablement entre les
            paires de classes, chacune apportant ses meilleures variables.

        Les doublons sont retires en conservant l'ordre de priorite.
        """
        k = int(min(k, X.shape[1]))
        scores_par_paire = self._scores_1d_par_paire(X)
        agrege = self._agreger(scores_par_paire)

        retenus = []
        vus = set()

        def ajouter(indices):
            for i in indices:
                i = int(i)
                if i not in vus:
                    vus.add(i)
                    retenus.append(i)

        part_agregee = max(1, k // 3)
        ajouter(np.argsort(agrege, kind="stable")[:part_agregee])
        quota = max(1, -(-(k - part_agregee) // len(self.paires)))
        for ligne in scores_par_paire:
            ajouter(np.argsort(ligne, kind="stable")[:quota])

        if len(retenus) < k:  # complement eventuel par le classement agrege
            ajouter(np.argsort(agrege, kind="stable")[:k])
        return retenus[:k]

    def evaluer_2d(self, X, j, k):
        valeurs = [overlap_2d(X, j, k, self.indices[c1], self.indices[c2])
                   for c1, c2 in self.paires]
        if any(not np.isfinite(v) for v in valeurs):
            return np.inf
        return float(self._agreger(np.array(valeurs)))

    def evaluer_colonne_fixee(self, X, k):
        return self._agreger(np.vstack([
            overlap_2d_colonne_fixee(X, k, self.indices[c1], self.indices[c2])
            for c1, c2 in self.paires
        ]))

    def detail(self, X, j, k):
        valeurs = {}
        for c1, c2 in self.paires:
            valeurs["{}-{}".format(c1, c2)] = overlap_2d(
                X, j, k, self.indices[c1], self.indices[c2])
        tableau = np.array(list(valeurs.values()), dtype=np.float64)
        finis = tableau[np.isfinite(tableau)]
        return {
            "overlap_max": float(tableau.max()),
            "overlap_moyen": float(finis.mean()) if finis.size else np.inf,
            "nb_paires_bien_separees": int((tableau < 1.0).sum()),
            "nb_paires": int(tableau.size),
            "par_paire": valeurs,
        }


# -----------------------------------------------------------------------------
# Strategie de recherche
# -----------------------------------------------------------------------------
def rechercher_meilleure_paire(X, evaluateur,
                               n_candidates=60,
                               max_iterations=10,
                               epsilon=1e-6,
                               exhaustif_si_p_max=20,
                               journal=None):
    """Recherche la meilleure paire de variables selon le schema de la p.37.

    Parametres
    ----------
    X : ndarray (n, p) - matrice des attributs candidats
    evaluateur : objet exposant `evaluer_1d`, `evaluer_2d` et
        `evaluer_colonne_fixee`
    n_candidates : taille de la liste de candidats apres filtrage univarie
    max_iterations : critere d'arret 1
    epsilon : amelioration relative minimale, critere d'arret 2
    exhaustif_si_p_max : si p est inferieur ou egal a cette valeur, une
        recherche exhaustive est effectuee (cas de l'ACP, 15 composantes).

    Retour
    ------
    dictionnaire decrivant la paire retenue et le deroulement de la recherche.
    """
    X = np.asarray(X, dtype=np.float64)
    p = X.shape[1]
    trace = []

    def noter(etape, message):
        trace.append("[{}] {}".format(etape, message))
        if journal:
            journal("    {} : {}".format(etape, message))

    # --- Cas particulier : recherche exhaustive possible ---------------------
    if p <= exhaustif_si_p_max:
        noter("Strategie", "recherche exhaustive ({} paires)".format(p * (p - 1) // 2))
        meilleur_score = np.inf
        meilleure_paire = (0, min(1, p - 1))
        evaluations = 0
        for j in range(p):
            for k in range(j + 1, p):
                score = evaluateur.evaluer_2d(X, j, k)
                evaluations += 1
                if score < meilleur_score:
                    meilleur_score = score
                    meilleure_paire = (j, k)
        noter("Critere d'arret", "exploration complete terminee")
        return _resultat(evaluateur, X, meilleure_paire, meilleur_score,
                         evaluations, 0, "exhaustive", trace)

    # --- Etape A : filtrage univarie (Attributs -> sous-ensemble) ------------
    scores_1d = evaluateur.evaluer_1d(X)
    evaluations = p
    candidats = evaluateur.candidats(X, n_candidates)
    k_candidats = len(candidats)
    noter("Strategie",
          "filtrage univarie sur {} variables, {} candidats retenus "
          "(meilleur score 1D = {:.4f})".format(p, k_candidats, scores_1d[candidats[0]]))

    # --- Etape B : balayage des candidats contre toutes les variables --------
    # Pour chaque variable candidate k, une passe vectorisee donne le score de
    # TOUTES les paires (j, k) avec j parcourant les p variables. En k_candidats
    # passes, on evalue donc k_candidats x p paires, soit une fraction
    # importante de l'espace de recherche pour un cout tres faible.
    meilleur_score = np.inf
    j_courant, k_courant = int(candidats[0]), int(candidats[min(1, k_candidats - 1)])
    for c in candidats:
        scores = evaluateur.evaluer_colonne_fixee(X, c)
        evaluations += p
        meilleur_j = int(np.argmin(scores))
        if scores[meilleur_j] < meilleur_score:
            meilleur_score = float(scores[meilleur_j])
            j_courant, k_courant = meilleur_j, int(c)
    noter("Evaluation",
          "{} passes vectorisees ({} paires evaluees) ; paire initiale "
          "({}, {}) de score {:.4f}"
          .format(k_candidats, k_candidats * p, j_courant, k_courant, meilleur_score))

    # --- Etape C : recherche locale par echange de coordonnees ---------------
    # A chaque demi-iteration, une variable est fixee et la meilleure partenaire
    # est cherchee parmi les p variables (et non parmi les seuls candidats), ce
    # qui permet de sortir du sous-ensemble issu du filtrage univarie.
    iteration = 0
    raison_arret = "nombre maximal d'iterations atteint"
    while iteration < max_iterations:
        iteration += 1
        score_debut = meilleur_score
        marge = epsilon * max(abs(meilleur_score), 1.0)

        # C.1 : la seconde variable est fixee, on cherche la meilleure premiere.
        scores = evaluateur.evaluer_colonne_fixee(X, k_courant)
        evaluations += p
        candidat_j = int(np.argmin(scores))
        if scores[candidat_j] < meilleur_score - marge:
            meilleur_score = float(scores[candidat_j])
            j_courant = candidat_j

        # C.2 : la premiere variable est fixee, on cherche la meilleure seconde.
        scores = evaluateur.evaluer_colonne_fixee(X, j_courant)
        evaluations += p
        candidat_k = int(np.argmin(scores))
        if scores[candidat_k] < meilleur_score - marge:
            meilleur_score = float(scores[candidat_k])
            k_courant = candidat_k

        amelioration = score_debut - meilleur_score
        reference = max(abs(score_debut), 1.0)
        noter("Iteration {}".format(iteration),
              "paire ({}, {}), score {:.4f}, amelioration {:.3e}"
              .format(j_courant, k_courant, meilleur_score, amelioration))

        if amelioration <= epsilon * reference:
            raison_arret = ("aucune amelioration relative superieure a "
                            "epsilon = {:.1e}".format(epsilon))
            break

    noter("Critere d'arret", raison_arret)
    return _resultat(evaluateur, X, (j_courant, k_courant), meilleur_score,
                     evaluations, iteration, "filtrage + recherche locale", trace)


def _resultat(evaluateur, X, paire, score, evaluations, iterations, strategie, trace):
    j, k = int(paire[0]), int(paire[1])
    resultat = {
        "variable_1": j,
        "variable_2": k,
        "score": float(score),
        "mesure": evaluateur.nom,
        "cible": evaluateur.cible,
        "strategie": strategie,
        "nb_paires_evaluees": int(evaluations),
        "nb_iterations": int(iterations),
        "trace": trace,
    }
    resultat.update(evaluateur.detail(X, j, k))
    return resultat
