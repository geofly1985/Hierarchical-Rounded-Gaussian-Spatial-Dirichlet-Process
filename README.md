# Hierarchical-Rounded-Gaussian-Spatial-Dirichlet-Process
This repo contains the code of HRGSDP approach proposed in "A Bayesian Nonparametric model for textural pattern heterogeneity" by Li, X et all. 

The rgsdp.cpp, dp.h, samplers.R and cal_stat.R will be called in the simulation scripts. Be sure to change the file path accordingly.

The folder with name 's = *' contains the 100 simulation scripts each, each simulation script will generate a csv file that contains the chi sqaure test statistic, number of clusters as well as mis-assignment rate.

figure_simu.R will generate the figure5,6 in the manuscripts for simulation results.

Note: In each single simulaton, we also performed convergence diagnosis for each parameter (Geweke's diagnostic), if any parameter doesn't pass converge check, we change the seed.
