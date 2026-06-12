# ia-bioreprod2026-06
Résultat Atelier Reproducibilité - IA bioscripting 2026 - SC


TP réalisé avec Perplexity academic, en mode computer.
Espace de travail de Perplexity : 
https://www.perplexity.ai/spaces/iabioscriting-2026-EhFNwHfDR82X24fqCQasTw

Premier chat (en bas) : session de départ, finie quand j'ai épuisé le compte gratuit
Deuxième chat (en haut) quand j'ai basculé sur le forfait Pro pour teminer le TP


# fig2aheatmap

`fig2aheatmap` est un dépôt R reproductible pour recréer une heatmap de type Figure 2A à partir du tableau supplémentaire publié dans Kelliher et al. 2016. La publication indique que la Figure 2A de *S. cerevisiae* contient 1246 gènes périodiques après filtrage Lomb-Scargle, ordonnés selon leur temps de pic d'expression. [file:22][file:23]

## Objectif

Le projet lit le fichier Excel supplémentaire, sélectionne les gènes avec `LS_cutoff == "Yes"`, ordonne les lignes selon `Figure2A_order_peaktime`, calcule un z-score par gène et génère à la fois une matrice CSV et une heatmap PDF. Les colonnes utilisées par cette logique sont bien présentes dans la table fournie. [file:23]

## Structure du dépôt

- `R/figure2a.R` : fonctions principales.
- `scripts/run_figure2a.R` : exécution en ligne de commande.
- `tests/testthat/` : tests unitaires.
- `inst/extdata/` : fichier Excel d'exemple.
- `.github/workflows/` : intégration continue.

## Mise en place locale

```r
install.packages(c("renv", "devtools", "roxygen2", "BiocManager"))
BiocManager::install("ComplexHeatmap")
renv::init()
renv::settings$snapshot.type("explicit")
renv::snapshot()
```

L'usage de `renv::snapshot()` permet d'enregistrer précisément l'état des dépendances dans `renv.lock`, ce qui est recommandé pour les workflows R reproductibles. [web:15][web:42][web:48]

## Exécution

```bash
Rscript scripts/run_figure2a.R
```

Les sorties sont écrites dans `outputs/`.

## Vérifications qualité

```r
devtools::document()
devtools::check()
testthat::test_dir("tests/testthat")
lintr::lint_package()
```

L'écosystème `usethis` et `r-lib/actions` est la voie recommandée pour préparer un dépôt R moderne avec documentation et CI GitHub Actions. [web:35][web:40][web:43]
