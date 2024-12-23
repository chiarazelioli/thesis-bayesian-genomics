# Master Thesis Project
## Bayesian Genomic Prediction: a Comparative Study
This repository contains scripts and data for my thesis, where we investigate the performance across several genetic architectures of two Bayesian methods for genomic prediction in wild populations.

We here provide little background. 
Phenotypic traits **y** are modelled via marker-based regression models **y = Xb + Zu + Wd +e**, mixed effects models accounting for fixed (**Xb**), random genetic (**Zu**) and random environmental (**Wd**) effects. **Z** is the Nxp matrix of the SNPs, which is high-dimensional and the number of markers (sequenced genetic information) p greatly exceeds the number of individuals N. Vector **u** contains the _SNP-effects_. _Genomic prediction_ aims at making inference on **u**, as **Zu** represents the _breeding value_,
a parameter that quantifies the additive effect of inherited genes on the phenotype of an individual.

In particular, in wild populations studies, the number of genotyped individuals N is even smaller, and more complex phenomena take place, like gene-gene and gene-environment interactions, together with higher variability within a population. Thus, we look for dimensionality reduction techiques. 

To pursue our goal, we carry out a simulation study where we reproduce ten genetic architectures of phenotypic traits, spanning from oligogenic to polygenic, using real SNPs, and then we simulate phenotypes according to the marker-based regression model. We use a reduced version of it, which accounts for random additive genetic effects only (**y = Zu + e**). The goal with genomic prediction is to estimate the _breeding value_.


The two methods we compared are:
1. BayesR, currently the state-of-the-art in this field [Link Text](http://dx.doi.org/10.3168/jds.2011-5019). BayesR assumes that most SNPs have negligible or a very small effect on the expression of a trait, while the remaining SNPs affect the trait to different degrees. This corresponds to a Gaussian mixture prior on the SNP-effects, and a Dirichlet prior on the probability vector.
This method deals with the high-dimensionality of genomic data thanks to the sparsity introduced by the prior.
BayesR uses MCMC to sample from the full posterior of **u**, and follows a two-step approach, only recently handled with one single run via 'hibayes' Rpackage.

3. BPCRR, namely Bayesian Principal Components Ridge Regression, is a newly published method [Link Text](https://doi.org/10.1101/2024.06.01.596874), which first reduces dimensionality via PCA of the SNP matrix, and then sets a Gaussian prior on the _PC-effects_, indeed assuming they are i.i.d. . Indeed, we retain k < N out of p columns containing majority of the information, with a reduced matrix **Z'**. Also, the vector of SNP-effects becomes **u'**, vector of PC-effects. BPCRR numerically estimates **u'** via INLA and is therefore much faster for large datasets.


## Data
The dataset, "SNP_100k.rds", contains 100'000 SNPs for 3032 genotyped individuals. 
Contact stefamu@ntnu.no to download it.

## Scripts

### Functions
- BPCRR.R: perform genomic prediction by BPCRR using a fixed prior on PC variances
- BPCRR_default.R: perform genomic prediction by BPCRR using a default prior on PC variances
- BPCRR_fix_effects.R: perform genomic prediction by BPCRR accounting for fixed effects **Xb** in model **y = Xb Zu + e**
- BayesR.R: perform genomic prediction by BayesR
- SVD.R: perform singular value decomposition of the SNP matrix to extract PCs.

### Genomic prediction (and possibly additive variance estimation) using BPCRR and BayesR
- demo_sim.R: demonstration script (credits to janne.c.h.aspheim@ntnu.no). I recommend to give a look to this before using any other script
- main.R: optimized script which runs several seeds for reliable predictions and creates .csv files to store the results so that we may directly use plot them for visual inspection
- main_all_architectures.R: version of 'main.R' which performs genomic prediction for all the architectures set at the beginning. May require long time to run.

### Robustness Assessment - BPCRR
- vary_nPCs.R: genomic prediction by BPCRR when choosing different numbers of PCs
- increase_prior_var.R: genomic prediction by BPCRR when manually varying the fixed prior variance
- ridge_on_PCs.R: perform genomic prediction following a frequentist ridge regression on the PCs used for BPCRR.

**NOTE on PC scaling**: In all the scripts, 'demo_sim.R', 'main.R', 'main_all_pis.R', 'main_accuracy_inspection.R' (UPDATE), when preparing the principal components, one can either scale PCs variances to 1 (standardize), by setting "scale = TRUE" in line 'XX <- scale(XX, center = TRUE, scale = FALSE)', or simply centering them keeping the line as it is, and later scaling them w.r.t. the variance of the first PC.

