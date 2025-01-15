# A Low-Rank Bayesian Approach for Geoadditive Modeling
Bryan Sumalinab, Oswaldo Gressani, Niel Hens and Christel Faes

## About this repository
This repository contains the R codes used to generate results for the paper "A Low-Rank Bayesian Approach for Geoadditive Modeling" by Bryan Sumalinab, Oswaldo Gressani, Niel Hens, and Christel Faes.

## Simulation
The R codes used in the simulation study are located in the **Simulation** folder. Inside this folder are subfolders for different models and scenarios: **NB** (Negative Binomial model assumptions for count data), **Poisson** (Poisson model assumptions for count data), and **Gaussian** (Gaussian data).  Each of these main folders (NB, Poisson, Gaussian) includes a **functions** folder with the necessary functions to run the simulation codes and subfolders named **s1**, **s2**, and **s3**, representing the two-dimensional functions considered in the simulation.

The main functionalities for each model type are implemented in the scripts ***Krig.NB*** (for Negative Binomial models), ***Krig.Pois*** (for Poisson models), and ***Krig.Gauss*** (for Gaussian models).  These scripts include an option *cov.model* that allows for the specification of different covariance functions, such as exponential, Matern, circular, and spherical.

The **Simulation** folder also contains simulation results for the INLA-SPDE approach (**INLA-SPDE**) and a comparison of the proposed low-rank Bayesian approach with the classical kriging approach (**Comparison with likfit**).

## Real Data Application
For real data applications, the **RealData** folder contains R codes to analyze the Meuse river data demonstrating the practical application of the proposed methods on real-world datasets.
