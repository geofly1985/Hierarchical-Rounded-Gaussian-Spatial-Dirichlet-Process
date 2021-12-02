# Hierarchical-Rounded-Gaussian-Spatial-Dirichlet-Process
This repo contains the code of HRGSDP approach proposed in "A Bayesian Nonparametric model for textural pattern heterogeneity" by Li, X et al. 
https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssc.12469?casa_token=gyBe-YjHP_sAAAAA:bwtdyHWUVYsOHTyvMryhz7q1FvSLlhVI7xRtnwvj8PmADPzJMGL5tF0Gbao59SF1q7QtrnCLsB-IDD4

The rgsdp.cpp, dp.h, samplers.R and cal_stat.R will be called in the simulation scripts. Be sure to change the file path accordingly.

The folder with name 's = *' contains the 100 simulation scripts each, each simulation script will (a) generate a dataset contains 100 GLCMs (with 5 patterns); (b) call the samplers to estimate the parameters; (c) The results will be stored in the log file which contains the chi sqaure test statistic, number of clusters as well as mis-assignment rate.

figure_simu.R will generate the figure 5, 6 in the manuscripts for simulation results.

Note: In each single simulaton, we also performed convergence diagnosis for each parameter (Geweke's diagnostic), if any parameter doesn't pass converge check, we change the seed.
