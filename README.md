# metabotyping
This repository contains the R and stan code used to perform the analysis discussed in the manuscript ["Metabotyping as a Stopover in Genome-to-Phenome Mapping" by Pubudu Handakumbura et al.](https://www.nature.com/articles/s41598-019-38483-0) published in Scientific Reports.  The data needed to run this code is available from osti.gov [here](https://www.osti.gov/biblio/1501917).

The contents of the `Code` folder are as follows:

  * bayes_allometric.R - This code fits the allometric model via Bayesian hierarchical modelling (BHM) using `stan`
  * BHM_allom_lscale.stan - This script defines the Bayesian hierarchical model
  * Create_figures_2_and_Ex2.R - Code necessary to create Figure 2 and Extended Data Figure 2.
  * Create_RF_Imp_Plot.R - This script collects the random forest results (created in 'rf_biomass_pred.R') and plots the features based on their feature importance
  * Dendrogram_Comparison.R - This script compares the above and belowground dendrograms using Baker's gamma then performs a permutation test to assess the statistical significance of the observed Baker's gamma value
  * metab_analysis.R - This script uses the `MSomicsSTAT` package to analyze the raw metabolite abundance data.  This script creates the heat maps in Figure 4, Extended Data Figure 4 and the scrollable heatmaps in [this shiny app](https://ascm.shinyapps.io/BAS_gobrachy/)
  * PCA_PERMONOVA.R - This script performs a principal component analysis and PERMONOMVA
  * pheno_funs.R - This script contains a variety of functions used in the other scripts
  * rf_biomass_pred.R - This script fits the random forest model to predict biomass as a function of metabolite abundances
  
