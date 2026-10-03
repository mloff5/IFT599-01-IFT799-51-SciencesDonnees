# IFT599 - Science de données - TP1

## Compréhension et visualisation des données (Automne 2026)

**Classe : IFT599-1**

**Équipe :**

| Nom | CIP |
|---|---|
| Mael Hounkou | houg1304 |
| Mouhamadou Moustapha Gueye | guem8961 |
| Cheikh Ahmadou Tidiane Ly | lyzc3821 |

---

## 1. Objet du programme

Le programme réalise l'étude exploratoire de séparabilité demandée dans l'énoncé
du TP1, sur les profils génomiques RNA-Seq (HiSeq) PANCAN de 801 patients
répartis en 5 types de tumeurs (BRCA, KIRC, COAD, LUAD, PRAD).

Il implante les deux méthodes exigées :

* **Méthode 1 — analyse de séparabilité sans visualisation.** Les mesures de
  cohésion (`dist_intra`), de séparation (`dist_inter`) et l'indicateur de
  superposition `Overlap()` sont calculées sur **toutes** les variables de
  l'ensemble de données de départ, avec la distance **euclidienne** et la
  distance de **Mahalanobis**, et présentées sous forme de tableaux.
* **Méthode 2 — analyse avec visualisation.** Sélection automatique de la
  meilleure paire de variables (variables originales et composantes
  principales) selon le schéma de sélection de dimensions de la diapositive 37
  du Thème 1, puis production des nuages de points avec distributions
  marginales. La projection t-SNE est également traitée.

---

## 2. Installation

Python 3.9 ou plus récent.

```bash
pip install -r requirements.txt
```

Dépendances : `numpy`, `pandas`, `matplotlib`, `scikit-learn`.

---

## 3. Mise en place des données

Les données **ne sont pas incluses** dans la remise (consigne de l'énoncé).
Décompressez l'archive `gene+expression+cancer+rna+seq.zip`, puis
`TCGA-PANCAN-HiSeq-801x20531.tar.gz`.

Le dossier `TCGA-PANCAN-HiSeq-801x20531/` (contenant `data.csv` et
`labels.csv`) est **localisé automatiquement** s'il se trouve à l'un des
emplacements suivants, relativement au dossier du code :

```
donnees/TCGA-PANCAN-HiSeq-801x20531/   <- projet autonome
TCGA-PANCAN-HiSeq-801x20531/           <- données à côté du code
../TCGA-PANCAN-HiSeq-801x20531/        <- code dans un sous-dossier
donnees/                               <- data.csv et labels.csv à plat
```

Les deux premiers couvrent le cas d'une utilisation autonome et celui d'une
intégration dans un dépôt où le jeu de données est déjà rangé à côté du code.

Il est toujours possible d'imposer un autre emplacement :

```bash
python main.py --donnees /chemin/vers/TCGA-PANCAN-HiSeq-801x20531
```

---

## 4. Exécution

```bash
python main.py                      # exécution complète (≈ 3 à 5 minutes)
python main.py --etape m1           # Méthode 1 seulement
python main.py --etape m2           # Méthode 2 seulement
python main.py --n-variables 1000   # autre taille d'ensemble de départ
python main.py --mode-selection variance   # sélection par variance au lieu d'aléatoire
python main.py --sans-histogrammes  # nuages de points sans marginales
python main.py --sans-cache         # force la relecture du fichier CSV
```

Le premier appel lit le fichier `data.csv` (≈ 200 Mo) et met l'ensemble de
départ en cache dans `cache/` ; les appels suivants sont immédiats.

---

## 5. Résultats produits

```
resultats/
├── SYNTHESE.md                      <- tous les tableaux, en Markdown
├── variables_selectionnees.txt      <- liste des dimensions retenues
├── tableaux/
│   ├── methode1_cohesion.csv
│   ├── methode1_separabilite.csv
│   ├── methode1_overlap_synthese.csv
│   ├── methode1_diagnostic_covariance.csv
│   ├── methode2_selection.csv
│   └── methode2_trace_selection.txt <- déroulement de la recherche (p.37)
└── figures/
    ├── methode1_overlap_<transformation>_<métrique>.png
    ├── acp_variance_expliquee.png
    ├── originales_brutes_*.png
    ├── originales_zscore_*.png
    ├── acp_*.png
    └── tsne_*.png
```

`SYNTHESE.md` est le point d'entrée : il rassemble l'ensemble des tableaux
chiffrés dans l'ordre où ils sont attendus dans le rapport.

---

## 6. Choix de conception

### 6.1 Réduction dimensionnelle initiale (section 2.2 de l'énoncé)

1. Retrait des variables **constantes** (écart-type nul) : 20 531 colonnes au
   départ, **20 264** variables non constantes. Ces variables ne peuvent
   contribuer à aucune mesure de distance et dégradent le conditionnement de la
   matrice de covariance.
2. **Tirage aléatoire** (graine 42) de **2 000** variables parmi les variables
   non constantes. L'énoncé autorise explicitement une sélection aléatoire.
   Le mode `--mode-selection variance` permet de comparer avec une sélection par
   variance décroissante.

La liste exacte des variables retenues est écrite dans
`resultats/variables_selectionnees.txt`, conformément à la note 3 de l'énoncé.

> Remarque : l'énoncé mentionne 16 382 variables dans l'ensemble initial. Le
> fichier `data.csv` distribué en contient 20 531, dont 20 264 non constantes.
> L'écart provient vraisemblablement d'un filtrage préalable différent ; il est
> sans incidence puisque la réduction est de toute façon à notre charge.

### 6.2 Distance de Mahalanobis

La distance de Mahalanobis est calculée par **blanchiment** : si
`Σ⁻¹ = Wᵀ W`, alors `d_M(x, y) = ‖W x − W y‖₂`. La matrice `W` est obtenue une
seule fois par décomposition spectrale, puis toutes les distances sont des
distances euclidiennes dans l'espace blanchi. Cela évite de recalculer une forme
quadratique pour chacun des couples de points.

Avec n = 801 patients et p = 2 000 variables, la covariance empirique est de
rang au plus 800 : elle n'est **pas inversible**. Une régularisation de type
ridge est appliquée (`Σ + λ·tr(Σ)/p·I`, λ = 10⁻²). La variante pseudo-inverse
est disponible via `config.METHODE_INVERSION_COVARIANCE = "pinv"`.

La covariance utilisée est la covariance **globale** du jeu de données : une
covariance propre à chaque classe rendrait les distances non comparables d'une
classe à l'autre.

### 6.3 Sélection automatique de variables (Méthode 2)

Le schéma de la diapositive 37 du Thème 1 est implanté intégralement :

```
Attributs → Stratégie de recherche → Sous-ensemble d'attributs
                  ↑                           ↓
            Critère d'arrêt  ←  ←  ←  ←  Évaluation
                  ↓ oui
         Attributs sélectionnés → Processus de validation
```

* **Stratégie de recherche** : filtrage univarié par `Overlap()` en dimension 1
  (p évaluations), puis recherche locale par échange de coordonnées
  (*hill-climbing*) sur les `N_CANDIDATES` meilleures variables. Le budget total
  reste de l'ordre de 2p évaluations, conformément à la note 2 de l'énoncé, au
  lieu des p(p−1)/2 évaluations d'une recherche exhaustive.
* **Évaluation** : deux mesures distinctes, toutes deux construites à partir de
  `Overlap()`, comme l'exige l'énoncé :
  * pour **une** paire de classes : `f(j,k) = Overlap(C_a, C_b)` ;
  * pour **toutes** les paires de classes : `g_max(j,k) = max Overlap` (critère
    minimax) et `g_moyenne(j,k) = moyenne des Overlap`. Les deux agrégations
    sont exécutées et comparées.
* **Critère d'arrêt** : aucune amélioration relative supérieure à ε = 10⁻⁶, ou
  nombre maximal d'itérations atteint.
* **Validation** : vérification de la condition `Overlap < 1` et production de
  la figure correspondante.

Pour l'ACP (15 composantes), p ≤ 20 et la recherche devient **exhaustive**
(105 paires), ce qui garantit l'optimalité dans cet espace.

### 6.4 Distance utilisée dans la Méthode 2

Distance **euclidienne**, comme conseillé dans l'énoncé pour des raisons de
simplicité.

### 6.5 Affichage 1D (point 2.b, optionnel)

Les nuages de points sont accompagnés des **histogrammes marginaux** de chaque
classe le long des deux axes, ce qui illustre la séparation des classes variable
par variable. Désactivable avec `--sans-histogrammes`.

---

## 7. Structure du code

```
TP1_IFT599/
├── main.py                 programme principal (enchaînement des étapes)
├── config.py               paramètres centralisés
├── requirements.txt
├── SOURCES.md              sources et références utilisées
├── README.md
└── src/
    ├── data_loader.py      chargement et réduction dimensionnelle initiale
    ├── distances.py        distance euclidienne et de Mahalanobis
    ├── separabilite.py     dist_intra, dist_inter, Overlap()
    ├── transformations.py  Z-score, min-max, ACP (maison), t-SNE
    ├── selection.py        sélection automatique de paires de variables
    ├── visualisation.py    nuages de points et histogrammes marginaux
    ├── methode1.py         Méthode 1 (tableaux)
    ├── methode2.py         Méthode 2 (sélection + figures)
    └── synthese.py         mise en forme des résultats
```

Tous les modules sont indépendants et documentés ; chaque fichier porte
l'identification de l'équipe et de la classe exigée par l'énoncé.

---

## 8. Sources

Voir `SOURCES.md`.
