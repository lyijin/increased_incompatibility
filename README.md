# Increased incompatibility of heterologous algal symbionts under thermal stress in the cnidarian-dinoflagellate model Aiptasia #

This repository documents code used in a manuscript that still under review.

## Description of folder contents ##

``differential_gene_expression/`` carries out differential gene expression across 3 strains x 2 temperatures (25 C, 32 C) x 4 replicates = 24 individual samples. DE was performed using `kallisto` and `sleuth`. Analyses were then carried out on the output tables with the R/Rmd scripts, producing **Fig. 1**. Code written by [@1cuigx](https://github.com/1cuigx).

``heatmap_expr_change/`` investigates the effect of heat stress on the expression of symbiosis-associated genes (Cui et al., 2019). The changes in expression was contrasted across CC7, CC7-B01 and H2 strains. As there were no Python code that could cleanly create what we envisioned, we plotted the upregulated genes and downregulated genes separately, and merged both halves with a vector graphics editing software (Illustrator). Material in this folder is the basis for **Fig. 4**. Code written by [@lyijin](https://github.com/lyijin).
