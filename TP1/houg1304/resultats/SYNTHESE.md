# IFT599 - Science de donnees - TP1
## Comprehension et visualisation des donnees (Automne 2026)

**Classe :** IFT599-1

**Equipe :**

- Mael Hounkou (houg1304)
- Mouhamadou Moustapha Gueye (guem8961)
- Cheikh Ahmadou Tidiane Ly (lyzc3821)

---

> Ce fichier rassemble tous les resultats numeriques produits par le
> programme. Les tableaux sont directement reutilisables dans le
> rapport. Les figures correspondantes se trouvent dans
> `resultats/figures/`.

## 1. Ensemble de donnees de depart

```
Patients (n)                      : 801
Variables du fichier brut         : 20531
Variables non constantes          : 20264
Variables de l'ensemble de depart : 2000
Classes                           : BRCA (300), COAD (78), KIRC (146), LUAD (141), PRAD (136)
```

La reduction de dimensionnalite initiale a consiste a retirer les variables constantes puis a tirer au hasard (graine 42) 2000 variables parmi les 20264 variables non constantes. La liste complete des variables retenues se trouve dans `resultats/variables_selectionnees.txt`.

## 2. Methode 1 - separabilite sans visualisation

Mesures calculees sur **toutes** les 2000 variables de l'ensemble de depart.

### 2.1 Synthese des valeurs de Overlap

Une valeur inferieure a 1 indique que les deux classes sont bien separees.

| paire | brut + euclidienne | brut + mahalanobis | minmax + euclidienne | minmax + mahalanobis | zscore + euclidienne | zscore + mahalanobis |
|---|---|---|---|---|---|---|
| BRCA / COAD | 1.3577 | 1.0320 | 1.5494 | 1.0176 | 1.7733 | 1.0163 |
| BRCA / KIRC | 1.2199 | 1.0506 | 1.4849 | 1.0343 | 1.8734 | 1.0318 |
| BRCA / LUAD | 1.4615 | 1.0390 | 1.6631 | 1.0235 | 1.9191 | 1.0214 |
| BRCA / PRAD | 1.3947 | 1.0606 | 1.5288 | 1.0317 | 1.9108 | 1.0280 |
| COAD / KIRC | 1.1730 | 1.0406 | 1.3958 | 1.0254 | 1.6722 | 1.0235 |
| COAD / LUAD | 1.3952 | 1.0230 | 1.4921 | 1.0124 | 1.6777 | 1.0118 |
| COAD / PRAD | 1.0626 | 1.0504 | 1.1905 | 1.0228 | 1.4881 | 1.0198 |
| KIRC / LUAD | 1.2486 | 1.0477 | 1.4222 | 1.0313 | 1.6997 | 1.0286 |
| KIRC / PRAD | 1.0982 | 1.0580 | 1.3512 | 1.0312 | 1.6984 | 1.0294 |
| LUAD / PRAD | 1.2746 | 1.0578 | 1.3920 | 1.0287 | 1.7259 | 1.0249 |

### 2.2 Cohesion (distance intra-classe)

| classe | brut + euclidienne | brut + mahalanobis | zscore + euclidienne | zscore + mahalanobis | minmax + euclidienne | minmax + mahalanobis |
|---|---|---|---|---|---|---|
| BRCA | 83.9508 | 28.0352 | 73.1450 | 28.1324 | 10.3205 | 28.0980 |
| COAD | 80.6111 | 27.7997 | 64.4168 | 27.9660 | 9.1665 | 27.9048 |
| KIRC | 83.4904 | 28.0100 | 77.9048 | 28.1046 | 10.5982 | 28.0705 |
| LUAD | 79.6258 | 28.0012 | 63.3840 | 28.0756 | 8.9802 | 28.0521 |
| PRAD | 71.1906 | 28.0595 | 71.1317 | 28.1234 | 8.5975 | 28.0517 |

### 2.3 Detail par combinaison

#### brut + distance euclidienne

| classe_1 | classe_2 | dist_intra_C1 | dist_intra_C2 | dist_inter | overlap | bien_separees |
|---|---|---|---|---|---|---|
| BRCA | COAD | 83.9508 | 80.6111 | 60.6033 | 1.3577 | False |
| BRCA | KIRC | 83.9508 | 83.4904 | 68.6308 | 1.2199 | False |
| BRCA | LUAD | 83.9508 | 79.6258 | 55.9608 | 1.4615 | False |
| BRCA | PRAD | 83.9508 | 71.1906 | 55.6188 | 1.3947 | False |
| COAD | KIRC | 80.6111 | 83.4904 | 69.9498 | 1.1730 | False |
| COAD | LUAD | 80.6111 | 79.6258 | 57.4257 | 1.3952 | False |
| COAD | PRAD | 80.6111 | 71.1906 | 71.4295 | 1.0626 | False |
| KIRC | LUAD | 83.4904 | 79.6258 | 65.3220 | 1.2486 | False |
| KIRC | PRAD | 83.4904 | 71.1906 | 70.4227 | 1.0982 | False |
| LUAD | PRAD | 79.6258 | 71.1906 | 59.1609 | 1.2746 | False |

Paires bien separees : **0/10**

#### brut + distance mahalanobis

| classe_1 | classe_2 | dist_intra_C1 | dist_intra_C2 | dist_inter | overlap | bien_separees |
|---|---|---|---|---|---|---|
| BRCA | COAD | 28.0352 | 27.7997 | 27.0524 | 1.0320 | False |
| BRCA | KIRC | 28.0352 | 28.0100 | 26.6733 | 1.0506 | False |
| BRCA | LUAD | 28.0352 | 28.0012 | 26.9666 | 1.0390 | False |
| BRCA | PRAD | 28.0352 | 28.0595 | 26.4451 | 1.0606 | False |
| COAD | KIRC | 27.7997 | 28.0100 | 26.8159 | 1.0406 | False |
| COAD | LUAD | 27.7997 | 28.0012 | 27.2735 | 1.0230 | False |
| COAD | PRAD | 27.7997 | 28.0595 | 26.5884 | 1.0504 | False |
| KIRC | LUAD | 28.0100 | 28.0012 | 26.7299 | 1.0477 | False |
| KIRC | PRAD | 28.0100 | 28.0595 | 26.4983 | 1.0580 | False |
| LUAD | PRAD | 28.0012 | 28.0595 | 26.4982 | 1.0578 | False |

Paires bien separees : **0/10**

#### zscore + distance euclidienne

| classe_1 | classe_2 | dist_intra_C1 | dist_intra_C2 | dist_inter | overlap | bien_separees |
|---|---|---|---|---|---|---|
| BRCA | COAD | 73.1450 | 64.4168 | 38.7865 | 1.7733 | False |
| BRCA | KIRC | 73.1450 | 77.9048 | 40.3144 | 1.8734 | False |
| BRCA | LUAD | 73.1450 | 63.3840 | 35.5720 | 1.9191 | False |
| BRCA | PRAD | 73.1450 | 71.1317 | 37.7536 | 1.9108 | False |
| COAD | KIRC | 64.4168 | 77.9048 | 42.5544 | 1.6722 | False |
| COAD | LUAD | 64.4168 | 63.3840 | 38.0887 | 1.6777 | False |
| COAD | PRAD | 64.4168 | 71.1317 | 45.5441 | 1.4881 | False |
| KIRC | LUAD | 77.9048 | 63.3840 | 41.5632 | 1.6997 | False |
| KIRC | PRAD | 77.9048 | 71.1317 | 43.8764 | 1.6984 | False |
| LUAD | PRAD | 63.3840 | 71.1317 | 38.9702 | 1.7259 | False |

Paires bien separees : **0/10**

#### zscore + distance mahalanobis

| classe_1 | classe_2 | dist_intra_C1 | dist_intra_C2 | dist_inter | overlap | bien_separees |
|---|---|---|---|---|---|---|
| BRCA | COAD | 28.1324 | 27.9660 | 27.5983 | 1.0163 | False |
| BRCA | KIRC | 28.1324 | 28.1046 | 27.2529 | 1.0318 | False |
| BRCA | LUAD | 28.1324 | 28.0756 | 27.5143 | 1.0214 | False |
| BRCA | PRAD | 28.1324 | 28.1234 | 27.3622 | 1.0280 | False |
| COAD | KIRC | 27.9660 | 28.1046 | 27.3917 | 1.0235 | False |
| COAD | LUAD | 27.9660 | 28.0756 | 27.6950 | 1.0118 | False |
| COAD | PRAD | 27.9660 | 28.1234 | 27.4991 | 1.0198 | False |
| KIRC | LUAD | 28.1046 | 28.0756 | 27.3078 | 1.0286 | False |
| KIRC | PRAD | 28.1046 | 28.1234 | 27.3116 | 1.0294 | False |
| LUAD | PRAD | 28.0756 | 28.1234 | 27.4163 | 1.0249 | False |

Paires bien separees : **0/10**

#### minmax + distance euclidienne

| classe_1 | classe_2 | dist_intra_C1 | dist_intra_C2 | dist_inter | overlap | bien_separees |
|---|---|---|---|---|---|---|
| BRCA | COAD | 10.3205 | 9.1665 | 6.2884 | 1.5494 | False |
| BRCA | KIRC | 10.3205 | 10.5982 | 7.0437 | 1.4849 | False |
| BRCA | LUAD | 10.3205 | 8.9802 | 5.8026 | 1.6631 | False |
| BRCA | PRAD | 10.3205 | 8.5975 | 6.1872 | 1.5288 | False |
| COAD | KIRC | 9.1665 | 10.5982 | 7.0801 | 1.3958 | False |
| COAD | LUAD | 9.1665 | 8.9802 | 6.0811 | 1.4921 | False |
| COAD | PRAD | 9.1665 | 8.5975 | 7.4607 | 1.1905 | False |
| KIRC | LUAD | 10.5982 | 8.9802 | 6.8832 | 1.4222 | False |
| KIRC | PRAD | 10.5982 | 8.5975 | 7.1033 | 1.3512 | False |
| LUAD | PRAD | 8.9802 | 8.5975 | 6.3137 | 1.3920 | False |

Paires bien separees : **0/10**

#### minmax + distance mahalanobis

| classe_1 | classe_2 | dist_intra_C1 | dist_intra_C2 | dist_inter | overlap | bien_separees |
|---|---|---|---|---|---|---|
| BRCA | COAD | 28.0980 | 27.9048 | 27.5161 | 1.0176 | False |
| BRCA | KIRC | 28.0980 | 28.0705 | 27.1538 | 1.0343 | False |
| BRCA | LUAD | 28.0980 | 28.0521 | 27.4318 | 1.0235 | False |
| BRCA | PRAD | 28.0980 | 28.0517 | 27.2135 | 1.0317 | False |
| COAD | KIRC | 27.9048 | 28.0705 | 27.2932 | 1.0254 | False |
| COAD | LUAD | 27.9048 | 28.0521 | 27.6357 | 1.0124 | False |
| COAD | PRAD | 27.9048 | 28.0517 | 27.3535 | 1.0228 | False |
| KIRC | LUAD | 28.0705 | 28.0521 | 27.2087 | 1.0313 | False |
| KIRC | PRAD | 28.0705 | 28.0517 | 27.2127 | 1.0312 | False |
| LUAD | PRAD | 28.0521 | 28.0517 | 27.2692 | 1.0287 | False |

Paires bien separees : **0/10**

### 2.4 Diagnostic de la matrice de covariance (distance de Mahalanobis)

| transformation | metrique | methode | lambda | delta_absolu | rang_numerique | dimension_blanchie | conditionnement |
|---|---|---|---|---|---|---|---|
| brut | mahalanobis | ridge | 0.010 | 0.018 | 800 | 2000 | 32666.612 |
| zscore | mahalanobis | ridge | 0.010 | 0.010 | 800 | 2000 | 21809.881 |
| minmax | mahalanobis | ridge | 0.010 | 0.000 | 800 | 2000 | 25921.057 |

Avec n = 801 patients et p = 2000 variables, la covariance empirique est de rang au plus 800 : elle est donc singuliere. La regularisation retenue est indiquee ci-dessus.

Remarque utile pour la discussion : la distance de Mahalanobis calculee avec la covariance exacte serait rigoureusement invariante a une transformation affine par variable (le Z-score et le min-max en font partie). Les ecarts observes entre les colonnes `brut`, `zscore` et `minmax` en Mahalanobis proviennent donc uniquement de la regularisation, qui n'est pas invariante a l'echelle des variables.

## 3. Methode 2 - selection de variables et visualisation

### 3.1 Meilleure paire de variables pour **toutes** les paires de classes

Deux mesures d'evaluation agregees sont comparees : le maximum des Overlap (critere minimax, le plus exigeant) et leur moyenne (qualite d'ensemble de l'affichage).

| modalite | agregation | variable_1 | variable_2 | overlap_max | overlap_moyen | paires_bien_separees | nb_paires_evaluees | nb_iterations |
|---|---|---|---|---|---|---|---|---|
| variables originales (brutes) | max | gene_5256 | gene_18235 | 8.3909 | 6.6558 | 0/10 | 126000 | 1 |
| variables originales (brutes) | moyenne | gene_7896 | gene_18235 | 9.9144 | 4.8449 | 0/10 | 126000 | 1 |
| variables originales (Z-score) | max | gene_18235 | gene_5254 | 8.0567 | 5.3811 | 0/10 | 130000 | 2 |
| variables originales (Z-score) | moyenne | gene_1315 | gene_18235 | 8.6686 | 5.0135 | 0/10 | 126000 | 1 |
| composantes principales (ACP) | max | CP3 | CP12 | 33.7415 | 13.5945 | 2/10 | 105 | 0 |
| composantes principales (ACP) | moyenne | CP1 | CP3 | 39.8599 | 10.0083 | 0/10 | 105 | 0 |
| t-SNE (projection globale) | max et moyenne | t-SNE 1 | t-SNE 2 | 2.5163 | 0.6323 | 8/10 | 0 | 0 |

### 3.2 Meilleure paire de variables pour **chaque** paire de classes

#### variables originales (brutes)

| paire de classes | variable_1 | variable_2 | Overlap | nb_paires_evaluees | nb_iterations |
|---|---|---|---|---|---|
| BRCA / COAD | gene_1866 | gene_7965 | 0.2826 | 126000 | 1 |
| BRCA / KIRC | gene_15473 | gene_7896 | 0.8916 | 126000 | 1 |
| BRCA / LUAD | gene_3339 | gene_7896 | 2.4043 | 126000 | 1 |
| BRCA / PRAD | gene_19913 | gene_9177 | 0.8754 | 126000 | 1 |
| COAD / KIRC | gene_14857 | gene_12013 | 0.4301 | 126000 | 1 |
| COAD / LUAD | gene_7965 | gene_12013 | 0.8829 | 126000 | 1 |
| COAD / PRAD | gene_3119 | gene_11910 | 0.3193 | 126000 | 1 |
| KIRC / LUAD | gene_12078 | gene_439 | 0.8183 | 126000 | 1 |
| KIRC / PRAD | gene_14130 | gene_11910 | 0.3176 | 126000 | 1 |
| LUAD / PRAD | gene_9177 | gene_11910 | 0.3779 | 126000 | 1 |

Paires bien separees : **9/10**

#### variables originales (Z-score)

| paire de classes | variable_1 | variable_2 | Overlap | nb_paires_evaluees | nb_iterations |
|---|---|---|---|---|---|
| BRCA / COAD | gene_4333 | gene_7965 | 0.2855 | 126000 | 1 |
| BRCA / KIRC | gene_15473 | gene_7896 | 0.9186 | 126000 | 1 |
| BRCA / LUAD | gene_15895 | gene_357 | 1.4492 | 130000 | 2 |
| BRCA / PRAD | gene_19913 | gene_9111 | 0.8602 | 126000 | 1 |
| COAD / KIRC | gene_15807 | gene_12013 | 0.4183 | 126000 | 1 |
| COAD / LUAD | gene_12013 | gene_15983 | 0.7216 | 126000 | 1 |
| COAD / PRAD | gene_12377 | gene_11910 | 0.3325 | 126000 | 1 |
| KIRC / LUAD | gene_12078 | gene_439 | 0.8899 | 126000 | 1 |
| KIRC / PRAD | gene_1749 | gene_11910 | 0.3212 | 126000 | 1 |
| LUAD / PRAD | gene_9177 | gene_11910 | 0.3522 | 126000 | 1 |

Paires bien separees : **9/10**

#### composantes principales (ACP)

| paire de classes | variable_1 | variable_2 | Overlap | nb_paires_evaluees | nb_iterations |
|---|---|---|---|---|---|
| BRCA / COAD | CP4 | CP6 | 1.7914 | 105 | 0 |
| BRCA / KIRC | CP1 | CP11 | 2.3976 | 105 | 0 |
| BRCA / LUAD | CP5 | CP7 | 5.1063 | 105 | 0 |
| BRCA / PRAD | CP3 | CP13 | 1.7304 | 105 | 0 |
| COAD / KIRC | CP1 | CP5 | 1.0371 | 105 | 0 |
| COAD / LUAD | CP4 | CP15 | 1.8175 | 105 | 0 |
| COAD / PRAD | CP3 | CP6 | 0.5685 | 105 | 0 |
| KIRC / LUAD | CP1 | CP5 | 1.6220 | 105 | 0 |
| KIRC / PRAD | CP1 | CP3 | 1.1211 | 105 | 0 |
| LUAD / PRAD | CP3 | CP7 | 0.9190 | 105 | 0 |

Paires bien separees : **2/10**

#### t-SNE (projection par paire)

| paire de classes | variable_1 | variable_2 | Overlap | nb_paires_evaluees | nb_iterations |
|---|---|---|---|---|---|
| BRCA / COAD | t-SNE 1 | t-SNE 2 | 0.3758 | 0 | 0 |
| BRCA / KIRC | t-SNE 1 | t-SNE 2 | 0.3563 | 0 | 0 |
| BRCA / LUAD | t-SNE 1 | t-SNE 2 | 2.1881 | 0 | 0 |
| BRCA / PRAD | t-SNE 1 | t-SNE 2 | 0.4865 | 0 | 0 |
| COAD / KIRC | t-SNE 1 | t-SNE 2 | 0.2521 | 0 | 0 |
| COAD / LUAD | t-SNE 1 | t-SNE 2 | 0.5645 | 0 | 0 |
| COAD / PRAD | t-SNE 1 | t-SNE 2 | 0.2277 | 0 | 0 |
| KIRC / LUAD | t-SNE 1 | t-SNE 2 | 0.3211 | 0 | 0 |
| KIRC / PRAD | t-SNE 1 | t-SNE 2 | 0.2759 | 0 | 0 |
| LUAD / PRAD | t-SNE 1 | t-SNE 2 | 0.3258 | 0 | 0 |

Paires bien separees : **9/10**

### 3.3 ACP - variance expliquee

| composante | valeur_propre | variance_expliquee_% | variance_cumulee_% |
|---|---|---|---|
| CP1 | 218.361 | 10.904 | 10.904 |
| CP2 | 167.338 | 8.356 | 19.261 |
| CP3 | 160.718 | 8.026 | 27.287 |
| CP4 | 99.948 | 4.991 | 32.278 |
| CP5 | 83.522 | 4.171 | 36.449 |
| CP6 | 56.775 | 2.835 | 39.284 |
| CP7 | 51.652 | 2.579 | 41.863 |
| CP8 | 42.340 | 2.114 | 43.978 |
| CP9 | 33.394 | 1.668 | 45.645 |
| CP10 | 24.752 | 1.236 | 46.881 |
| CP11 | 20.482 | 1.023 | 47.904 |
| CP12 | 19.652 | 0.981 | 48.886 |
| CP13 | 17.842 | 0.891 | 49.777 |
| CP14 | 14.862 | 0.742 | 50.519 |
| CP15 | 13.934 | 0.696 | 51.215 |

## 4. Elements de comparaison attendus dans le rapport

Les tableaux ci-dessus permettent de traiter les trois comparaisons exigees au point 3 de la section 2.2 de l'enonce :

1. **Transformation vs non-transformation** : comparer les colonnes `brut`, `zscore` et `minmax` du tableau 2.1, ainsi que les modalites `variables originales (brutes)` et `variables originales (Z-score)` de la section 3.
2. **Euclidienne vs Mahalanobis** : comparer, a transformation egale, les colonnes `+ euclidienne` et `+ mahalanobis` du tableau 2.1, en tenant compte du diagnostic de la section 2.4.
3. **ACP vs t-SNE** : comparer les lignes correspondantes des tableaux 3.1 et 3.2 et les figures `acp_*.png` et `tsne_*.png`.

---

Fichier genere automatiquement par `main.py`.
