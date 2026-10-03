# Sources et références

**IFT599 - Science de données - TP1 (Automne 2026) - Classe : IFT599-1**

| Nom | CIP |
|---|---|
| Mael Hounkou | houg1304 |
| Mouhamadou Moustapha Gueye | guem8961 |
| Cheikh Ahmadou Tidiane Ly | lyzc3821 |

---

## 1. Jeu de données

* **RNA-Seq (HiSeq) PANCAN**, UCI Machine Learning Repository.
  <https://archive.ics.uci.edu/dataset/401/gene+expression+cancer+rna+seq>
  Extraction aléatoire d'expressions génétiques de 801 patients répartis en
  5 types de tumeurs (BRCA, KIRC, COAD, LUAD, PRAD).
  Également fourni dans Teams du cours (General/travaux pratiques/TP1).

## 2. Matériel de cours

* **Thème 1 — Concepts de base**, Shengrui Wang, IFT599/IFT799, Automne 2026.
  * Diapositive 34 : normalisation min-max.
  * Diapositive 35 : standardisation Z-score.
  * Diapositives 36-37 : sélection de dimensions, schéma général de la
    procédure de sélection (attributs → stratégie de recherche →
    sous-ensemble → évaluation → critère d'arrêt → attributs sélectionnés).
    Ce schéma est implanté dans `src/selection.py`.
  * Diapositives 38-39 : analyse en composantes principales (ACP). La procédure
    y est décrite étape par étape (centrage-réduction, matrice de
    variance-covariance ou de corrélation, valeurs et vecteurs propres,
    projection sur les axes principaux) et est reproduite dans la classe `ACP`
    de `src/transformations.py`.
  * Diapositive 40 : approches non linéaires de réduction de dimension (t-SNE).
* **Énoncé du TP1**, *TP1 : Compréhension et visualisation des données*,
  IFT599/IFT799, Automne 2026. Toutes les formules de cohésion
  (`dist_intra`), de séparation (`dist_inter`) et de l'indicateur de
  superposition `Overlap()` implantées dans `src/separabilite.py` proviennent
  de la section 1.2 de ce document.

## 3. Bibliothèques utilisées

| Bibliothèque | Usage dans le projet | Référence |
|---|---|---|
| **NumPy** | algèbre linéaire, décomposition spectrale, SVD, calculs vectorisés | Harris, C.R. et al. (2020). *Array programming with NumPy*. Nature 585:357-362. <https://numpy.org/> |
| **pandas** | lecture des fichiers CSV, construction et export des tableaux | McKinney, W. (2010). *Data Structures for Statistical Computing in Python*. <https://pandas.pydata.org/> |
| **Matplotlib** | production de toutes les figures | Hunter, J.D. (2007). *Matplotlib: A 2D Graphics Environment*. Computing in Science & Engineering 9(3):90-95. <https://matplotlib.org/> |
| **scikit-learn** | **uniquement** pour l'algorithme t-SNE (`sklearn.manifold.TSNE`) | Pedregosa, F. et al. (2011). *Scikit-learn: Machine Learning in Python*. JMLR 12:2825-2830. <https://scikit-learn.org/> |

## 4. Algorithme t-SNE

* van der Maaten, L. et Hinton, G. (2008). *Visualizing Data using t-SNE*.
  Journal of Machine Learning Research, 9:2579-2605.
* Documentation de l'implantation utilisée :
  <https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html>

Conformément à la section 2.1 de l'énoncé, l'algorithme t-SNE n'a pas été
reprogrammé : l'implantation de scikit-learn est utilisée et citée ici.

## 5. Mise en page des figures

* La présentation « nuage de points central avec histogrammes marginaux » est
  celle indiquée en exemple dans l'énoncé :
  <https://seaborn.pydata.org/generated/seaborn.jointplot.html>
  Elle a été reproduite directement avec `matplotlib.gridspec` dans
  `src/visualisation.py` ; aucun code de seaborn n'a été copié et la
  bibliothèque n'est pas une dépendance du projet.

## 6. Distance de Mahalanobis

La formule `d_M(x, y) = sqrt((x − y)ᵀ Σ⁻¹ (x − y))` est la définition standard
de la distance de Mahalanobis, mentionnée dans l'énoncé du TP1 (section 1.2,
point 2). Son implantation par blanchiment (`Σ⁻¹ = Wᵀ W` donc
`d_M(x, y) = ‖W x − W y‖₂`) ainsi que la régularisation de type ridge employée
pour contourner la singularité de la covariance (n = 801 < p = 2000) ont été
écrites par l'équipe dans `src/distances.py`.

## 7. Code écrit par l'équipe

Les éléments suivants ne proviennent d'aucune source externe et ont été
entièrement programmés par l'équipe :

* les mesures `dist_intra`, `dist_inter` et `Overlap()` (`src/separabilite.py`),
  ainsi que leur version vectorisée pour l'évaluation univariée ;
* l'ACP (`src/transformations.py`), implantée avec NumPy en suivant la
  procédure des diapositives 38-39 du Thème 1 ;
* la distance de Mahalanobis par blanchiment et sa régularisation
  (`src/distances.py`) ;
* l'algorithme de sélection automatique de paires de variables
  (`src/selection.py`) : filtrage univarié, recherche locale par échange de
  coordonnées, mesures d'évaluation agrégées et critères d'arrêt ;
* l'ensemble de l'orchestration, des tableaux et des figures.

Conformément à la section 2.1 de l'énoncé, aucune ressource Internet n'a été
utilisée pour la partie analyse des résultats.
