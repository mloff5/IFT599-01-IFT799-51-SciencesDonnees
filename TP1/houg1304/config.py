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
# Fichier : config.py
# Role    : parametres globaux du projet (chemins, graine aleatoire, options).
# =============================================================================
"""Parametres centralises du TP1.

Tout parametre modifiable par l'utilisateur est regroupe ici afin de ne jamais
avoir a toucher au code des modules de calcul.
"""

from pathlib import Path

# -----------------------------------------------------------------------------
# Identification de l'equipe (reprise dans les en-tetes de tous les resultats)
# -----------------------------------------------------------------------------
CLASSE = "IFT599-1"
EQUIPE = [
    ("Mael Hounkou", "houg1304"),
    ("Mouhamadou Moustapha Gueye", "guem8961"),
    ("Cheikh Ahmadou Tidiane Ly", "lyzc3821"),
]

# -----------------------------------------------------------------------------
# Chemins
# -----------------------------------------------------------------------------
RACINE = Path(__file__).resolve().parent

FICHIER_DATA = "data.csv"
FICHIER_LABELS = "labels.csv"

# Repertoire contenant data.csv et labels.csv (jeu TCGA-PANCAN-HiSeq-801x20531).
# Il est localise automatiquement parmi les emplacements usuels, ce qui permet
# d'utiliser le projet aussi bien seul (dossier "donnees/") qu'integre a un
# depot ou le jeu de donnees est deja range a cote du code. L'emplacement peut
# toujours etre impose en ligne de commande avec --donnees.
NOM_DOSSIER_DONNEES = "TCGA-PANCAN-HiSeq-801x20531"

EMPLACEMENTS_DONNEES = [
    RACINE / "donnees" / NOM_DOSSIER_DONNEES,   # projet autonome
    RACINE / NOM_DOSSIER_DONNEES,               # donnees a cote du code
    RACINE.parent / NOM_DOSSIER_DONNEES,        # code dans un sous-dossier
    RACINE / "donnees",                         # donnees a plat dans donnees/
    RACINE.parent / "donnees" / NOM_DOSSIER_DONNEES,
]


def localiser_donnees():
    """Retourne le premier emplacement contenant data.csv et labels.csv.

    Si aucun candidat ne convient, le premier emplacement de la liste est
    retourne : le message d'erreur du chargement indiquera alors ou placer les
    fichiers.
    """
    for candidat in EMPLACEMENTS_DONNEES:
        if (candidat / FICHIER_DATA).exists() and (candidat / FICHIER_LABELS).exists():
            return candidat
    return EMPLACEMENTS_DONNEES[0]


DOSSIER_DONNEES = localiser_donnees()

DOSSIER_RESULTATS = RACINE / "resultats"
DOSSIER_TABLEAUX = DOSSIER_RESULTATS / "tableaux"
DOSSIER_FIGURES = DOSSIER_RESULTATS / "figures"
DOSSIER_CACHE = RACINE / "cache"

# -----------------------------------------------------------------------------
# Reproductibilite
# -----------------------------------------------------------------------------
GRAINE = 42

# -----------------------------------------------------------------------------
# Reduction dimensionnelle initiale (section 2.2 de l'enonce)
# -----------------------------------------------------------------------------
# Les 20531 colonnes du fichier brut contiennent des genes constants (variance
# nulle) qui n'apportent aucune information : ils sont retires d'office.
# On tire ensuite au hasard N_VARIABLES genes parmi les genes non constants pour
# constituer l'ensemble de donnees de depart.
N_VARIABLES = 2000

# "aleatoire"  : tirage uniforme parmi les variables non constantes (conforme a
#                l'enonce, qui autorise explicitement une selection aleatoire).
# "variance"   : conserve les N_VARIABLES de plus forte variance (mode de
#                comparaison, utilise pour la discussion du rapport).
MODE_SELECTION_INITIALE = "aleatoire"

# -----------------------------------------------------------------------------
# Distance de Mahalanobis
# -----------------------------------------------------------------------------
# La matrice de covariance estimee sur 801 patients pour 2000 variables est
# necessairement singuliere (rang <= 800). Deux strategies sont disponibles :
#   "ridge"  : Sigma + lambda * (trace(Sigma)/p) * I  -> inversible, stable.
#   "pinv"   : pseudo-inverse de Moore-Penrose (blanchiment du sous-espace de
#              rang plein uniquement).
METHODE_INVERSION_COVARIANCE = "ridge"
LAMBDA_RIDGE = 1e-2

# -----------------------------------------------------------------------------
# Methode 1 : combinaisons transformation x metrique a evaluer
# -----------------------------------------------------------------------------
TRANSFORMATIONS_METHODE1 = ["brut", "zscore", "minmax"]
METRIQUES_METHODE1 = ["euclidienne", "mahalanobis"]

# -----------------------------------------------------------------------------
# Methode 2 : selection automatique de paires de variables
# -----------------------------------------------------------------------------
# Nombre de variables candidates conservees apres le filtrage univarie.
# Chaque candidat declenche une passe vectorisee evaluant ses p paires, soit
# N_CANDIDATES x p paires explorees (~6 % de la recherche exhaustive pour
# p = 2000), bien loin des p(p-1)/2 paires d'une exploration complete.
N_CANDIDATES = 60
# Nombre maximal d'iterations de la recherche par echange (critere d'arret 1).
MAX_ITERATIONS = 10
# Amelioration relative minimale pour poursuivre la recherche (critere d'arret 2).
EPSILON_ARRET = 1e-6

# Nombre de composantes principales generees (l'enonce exige au moins 15).
N_COMPOSANTES_ACP = 15

# Parametres t-SNE
TSNE_PERPLEXITE = 30.0
TSNE_N_ITER = 1000
TSNE_LEARNING_RATE = "auto"

# -----------------------------------------------------------------------------
# Figures
# -----------------------------------------------------------------------------
# Palette fixe par classe (identique dans toutes les figures du rapport).
COULEURS_CLASSES = {
    "BRCA": "#1f77b4",
    "KIRC": "#d62728",
    "COAD": "#2ca02c",
    "LUAD": "#ff7f0e",
    "PRAD": "#9467bd",
}
MARQUEURS_CLASSES = {
    "BRCA": "o",
    "KIRC": "s",
    "COAD": "^",
    "LUAD": "D",
    "PRAD": "v",
}
DPI_FIGURES = 150
# Affichage conjoint des histogrammes marginaux (point 2.b optionnel de l'enonce)
AFFICHER_HISTOGRAMMES = True


def entete_texte(titre: str) -> str:
    """En-tete commun place en tete de chaque fichier de resultats."""
    membres = "\n".join(
        "#   {:<32s} ({})".format(nom, cip) for nom, cip in EQUIPE
    )
    return (
        "# =============================================================\n"
        "# IFT599 - Science de donnees - TP1 (Automne 2026)\n"
        "# Classe : {classe}\n"
        "# Equipe :\n{membres}\n"
        "#\n"
        "# {titre}\n"
        "# =============================================================\n"
    ).format(classe=CLASSE, membres=membres, titre=titre)
