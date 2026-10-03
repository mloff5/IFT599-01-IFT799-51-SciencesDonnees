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
# Fichier : src/data_loader.py
# Role    : chargement du jeu de donnees RNA-Seq (HiSeq) PANCAN et reduction
#           dimensionnelle initiale demandee a la section 2.2 de l'enonce.
#
# Sources : pandas (lecture CSV), numpy (algebre lineaire).
#           Aucun code externe n'a ete copie dans ce module.
# =============================================================================
"""Chargement des donnees et construction de l'ensemble de depart.

Le fichier `data.csv` contient 801 patients x 20531 genes. L'enonce demande de
reduire la dimensionnalite des le depart. La procedure retenue est la suivante :

1. Retrait des genes constants (variance nulle sur les 801 patients) : ils ne
   peuvent contribuer a aucune mesure de separation et rendent la matrice de
   covariance degeneree.
2. Tirage aleatoire (graine fixe) de `N_VARIABLES` genes parmi les genes non
   constants. L'enonce autorise explicitement une selection aleatoire.

Le resultat est mis en cache (format .npz) afin que les executions suivantes
soient immediates.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import config  # noqa: E402


class JeuDeDonnees:
    """Conteneur simple pour l'ensemble de donnees de depart."""

    def __init__(self, X, y, noms_variables, noms_patients, p_initial, p_non_constant):
        self.X = X                            # (n, p) float64
        self.y = y                            # (n,) str
        self.noms_variables = noms_variables  # (p,) str
        self.noms_patients = noms_patients    # (n,) str
        self.p_initial = p_initial            # nb de colonnes du fichier brut
        self.p_non_constant = p_non_constant  # nb de variables non constantes

    @property
    def n(self):
        return self.X.shape[0]

    @property
    def p(self):
        return self.X.shape[1]

    @property
    def classes(self):
        """Classes triees par effectif decroissant (ordre stable)."""
        valeurs, effectifs = np.unique(self.y, return_counts=True)
        ordre = np.argsort(-effectifs)
        return [str(v) for v in valeurs[ordre]]

    def effectifs(self):
        valeurs, effectifs = np.unique(self.y, return_counts=True)
        return {str(v): int(c) for v, c in zip(valeurs, effectifs)}

    def indices_classe(self, classe):
        return np.flatnonzero(self.y == classe)

    def sous_ensemble(self, classe):
        return self.X[self.indices_classe(classe)]

    def paires_de_classes(self):
        """Toutes les paires mutuelles de classes, sans repetition."""
        cls = sorted(self.classes)
        return [(cls[i], cls[j])
                for i in range(len(cls))
                for j in range(i + 1, len(cls))]

    def resume(self):
        lignes = [
            "Patients (n)                      : {}".format(self.n),
            "Variables du fichier brut         : {}".format(self.p_initial),
            "Variables non constantes          : {}".format(self.p_non_constant),
            "Variables de l'ensemble de depart : {}".format(self.p),
            "Classes                           : {}".format(
                ", ".join("{} ({})".format(c, n)
                          for c, n in sorted(self.effectifs().items()))),
        ]
        return "\n".join(lignes)


def _chemin_cache(n_variables, graine, mode):
    nom = "ensemble_depart_p{}_g{}_{}.npz".format(n_variables, graine, mode)
    return config.DOSSIER_CACHE / nom


def charger(dossier_donnees=None,
            n_variables=None,
            graine=None,
            mode=None,
            utiliser_cache=True,
            verbeux=True) -> JeuDeDonnees:
    """Charge le jeu de donnees et construit l'ensemble de depart reduit."""
    dossier_donnees = Path(dossier_donnees or config.DOSSIER_DONNEES)
    n_variables = n_variables or config.N_VARIABLES
    graine = config.GRAINE if graine is None else graine
    mode = mode or config.MODE_SELECTION_INITIALE

    cache = _chemin_cache(n_variables, graine, mode)
    if utiliser_cache and cache.exists():
        if verbeux:
            print("[donnees] lecture du cache : {}".format(cache.name))
        with np.load(cache, allow_pickle=False) as z:
            return JeuDeDonnees(
                X=z["X"].astype(np.float64),
                y=z["y"].astype(str),
                noms_variables=z["noms_variables"].astype(str),
                noms_patients=z["noms_patients"].astype(str),
                p_initial=int(z["p_initial"]),
                p_non_constant=int(z["p_non_constant"]),
            )

    chemin_data = dossier_donnees / config.FICHIER_DATA
    chemin_labels = dossier_donnees / config.FICHIER_LABELS
    if not chemin_data.exists() or not chemin_labels.exists():
        emplacements = "\n".join(
            "  - {}".format(c) for c in config.EMPLACEMENTS_DONNEES)
        raise FileNotFoundError(
            "data.csv et labels.csv sont introuvables dans {}.\n"
            "Emplacements explores automatiquement :\n{}\n"
            "Placez le dossier {} dans l'un d'eux, ou indiquez son chemin avec "
            "l'option --donnees."
            .format(dossier_donnees, emplacements, config.NOM_DOSSIER_DONNEES)
        )

    if verbeux:
        print("[donnees] lecture de {} (~200 Mo, patientez)...".format(config.FICHIER_DATA))
    donnees = pd.read_csv(chemin_data, index_col=0)
    etiquettes = pd.read_csv(chemin_labels, index_col=0)

    # Alignement strict patients <-> etiquettes.
    etiquettes = etiquettes.loc[donnees.index]
    y = etiquettes["Class"].to_numpy().astype(str)
    noms_patients = donnees.index.to_numpy().astype(str)
    noms_colonnes = donnees.columns.to_numpy().astype(str)
    X_complet = donnees.to_numpy(dtype=np.float64)
    p_initial = X_complet.shape[1]

    # --- Etape 1 : retrait des variables constantes --------------------------
    ecarts = X_complet.std(axis=0)
    non_constantes = np.flatnonzero(ecarts > 0.0)
    p_non_constant = int(non_constantes.size)
    if verbeux:
        print("[donnees] {} variables sur {} sont non constantes."
              .format(p_non_constant, p_initial))

    # --- Etape 2 : reduction a n_variables -----------------------------------
    if n_variables >= p_non_constant:
        retenues = non_constantes
    elif mode == "variance":
        variances = X_complet[:, non_constantes].var(axis=0)
        meilleures = np.argsort(-variances)[:n_variables]
        retenues = np.sort(non_constantes[meilleures])
    elif mode == "aleatoire":
        rng = np.random.default_rng(graine)
        tirage = rng.choice(p_non_constant, size=n_variables, replace=False)
        retenues = np.sort(non_constantes[tirage])
    else:
        raise ValueError("Mode de selection initiale inconnu : {}".format(mode))

    X = np.ascontiguousarray(X_complet[:, retenues])
    noms_variables = noms_colonnes[retenues]

    jeu = JeuDeDonnees(X, y, noms_variables, noms_patients,
                       p_initial, p_non_constant)

    config.DOSSIER_CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache,
        X=X, y=y, noms_variables=noms_variables, noms_patients=noms_patients,
        p_initial=p_initial, p_non_constant=p_non_constant,
    )
    if verbeux:
        print("[donnees] ensemble de depart mis en cache : {}".format(cache.name))
    return jeu


def ecrire_liste_variables(jeu: JeuDeDonnees, chemin):
    """Ecrit la liste des dimensions choisies (exige par la note 3 de l'enonce)."""
    chemin = Path(chemin)
    chemin.parent.mkdir(parents=True, exist_ok=True)
    with open(chemin, "w", encoding="utf-8") as f:
        f.write(config.entete_texte(
            "Liste des {} variables (genes) de l'ensemble de donnees de depart"
            .format(jeu.p)))
        f.write("# Selection : mode '{}', graine {}\n".format(
            config.MODE_SELECTION_INITIALE, config.GRAINE))
        f.write("# Variables non constantes disponibles : {}\n\n".format(
            jeu.p_non_constant))
        for nom in jeu.noms_variables:
            f.write(nom + "\n")
    return chemin
