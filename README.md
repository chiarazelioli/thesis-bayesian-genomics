# Master Thesis Project
## Bayesian Genomic Prediction: a Comparative Study
This repository contains scripts and data for my thesis, where we investigate the performance across several genetic architectures of two Bayesian methods for genomic prediction in wild populations.

We here provide a bit of background. 
Phenotypic traits **y** are modelled via marker-based regression models **y = Xb + Zu + Wd +e**, mixed effects models accounting for fixed (**Xb**), random genetic (**Zu**) and random environmental (**Wd**) effects. **Z** is the Nxp matrix of the SNPs, which is high-dimensional and the number of markers (sequenced genetic information) p greatly exceeds the number of individuals N. Vector **u** contains the _SNP-effects_. _Genomic prediction_ aims at making inference on **u**, as **Zu** represents the _breeding value_,
a parameter that quantifies the additive effect of inherited genes on the phenotype of an individual.

In particular, in wild populations studies, the number of genotyped individuals N is even smaller, and more complex phenomena take place, like gene-gene and gene-environment interactions, together with higher variability within a population. Thus, we look for dimensionality reduction techiques. 

To pursue our goal, we carry out a simulation study where we reproduce ten genetic architectures of phenotypic traits, spanning from oligogenic to polygenic, using real SNPs, and then we simulate phenotypes according to the marker-based regression model. We use a reduced version of it, which accounts for random additive genetic effects only (**y = Zu + e**). The goal with genomic prediction is to estimate the _breeding value_.


The two methods we compared are:
1. BayesR, currently the state-of-the-art in this field [Link Text](http://dx.doi.org/10.3168/jds.2011-5019). BayesR assumes that most SNPs have negligible or a very small effect on the expression of a trait, while the remaining SNPs affect the trait to different degrees. This corresponds to a Gaussian mixture prior on the SNP-effects, and a Dirichlet prior on the probability vector.
This method deals with the high-dimensionality of genomic data thanks to the sparsity introduced by the prior.
BayesR uses MCMC to sample from the full posterior of **u**, and follows a two-step approach, only recently handled with one single run via 'hibayes' Rpackage.

2. BPCRR, namely Bayesian Principal Components Ridge Regression, is a newly published method [Link Text](https://doi.org/10.1101/2024.06.01.596874), which first reduces dimensionality via PCA of the SNP matrix, and then sets a Gaussian prior on the _PC-effects_, indeed assuming they are i.i.d. . Indeed, we retain k < N out of p columns containing majority of the information, with a reduced matrix **Z'**. Also, the vector of SNP-effects becomes **u'**, vector of PC-effects. BPCRR numerically estimates **u'** via INLA and is therefore much faster for large datasets.

We also include scripts we used to assess robustness of BPCRR.

## Data
The data, "SNP_100k.rds", contains 100'000 SNPs for 3032 genotyped individuals. 
Contact stefamu@ntnu.no to download it.

Remark: we used these data to simulate phenotypes, which are the responses and which we call simulated _datasets_.

## Scripts

### A. External Functions
We write outer functions for readability.
- BPCRR.R: perform genomic prediction by BPCRR using a fixed prior on PC variances.
- BPCRR_default.R: perform genomic prediction by BPCRR using a default prior on PC variances.
- BPCRR_fix_effects.R: perform genomic prediction by BPCRR accounting for fixed effects **Xb** in model **y = Xb + Zu + Wd + e**.
- BayesR.R: perform genomic prediction by BayesR.
- SVD.R: perform singular value decomposition of the SNP matrix to extract PCs.


### B. Script for Genomic Prediction using BPCRR and BayesR
Here we pursue the main goal, which is to predict anc compare breeding values given by BPCRR and BayesR, having a simulated phenotype as response variable, with a genetic architecture decided a priori.

One can use any of the following three scripts to do so.

- demo_sim.R: demonstration script (credits to janne.c.h.aspheim@ntnu.no) perfect to first understand the workflow. A user sets each time a seed for reproducibility and a genetic architecture.
- main.R: optimized script which runs five seeds to simulate five datasets for reliable predictions and creates '.csv' files to store the results so that we may directly plot them for visual inspection.
- main_all_architectures.R: generalization of 'main.R', performing genomic prediction for all the architectures set at the beginning, simulating five datasets each. May require long time to run. At the end of this script there is an example of plot.

I suggest to start with 'demo_sim.R' for some first attempts. Then rely on 'main.R'. Use 'main_all_architectures.R' only if you want to handle more architectures at once.

#### B1. Workflow overview:
1) Parameters Setting: User defines the seed for reproducibility and the genetic architecture probabilities
2) Data Preparation: Load and scale the SNP matrix (SNP)
3) Phenotypes Simulation and Cross-Validation set-up
4) Dimensional Reduction (SVD)
5) Genomic Prediction via BPCRR and BayesR

#### B2. Phenotypes Simulation details:
0) Genetic architecture is specified by probabilities π1,π2,π3,π4 of the SNPs to have an effect drawn from a zero-mean Gaussian ditribution with variance, respectively, 0, 10^(-4), 10^(-3), 10^(-2). For example, a highly oligogenic trait will have (π1,π2,π3,π4)=(0.99,0,0,0.01) to represent that 99% of the SNPs does not affect the expression of the phenotype **y**.
1) Assign Effect Sizes to SNPs:
   Each SNP is assigned an effect size based on the probabilities and corresponding variances, then effect sizes are drawn from the specified normal distributions and 
breeding values are generated as a linear combination of the SNP matrix and effect sizes (**Zu**). Breeding values are scaled so their variance equals 1.
2) Add Environmental Effects (**e**):
   Environmental effects are simulated as random noise from a normal distribution with variance 2.
3) Combine Breeding Values and Environmental Effects:
   We eventually obtain pnenotypes according to model **y = Zu + e**. This yelds that the variance of the breeding value (additive genetic variance) is approximately 0.33, while the residual variance contributes the rest, providing an upper bound of sqrt(0.33) for prediction accuracy, since the latter is computed as cor(simulated phenotype, predicted breeding value) in 'BPCRR.R' and 'BayesR.R'. 


### C. Scripts for Robustness Assessment - BPCRR
These scripts will be public soon:
- vary_nPCs.R: genomic prediction by BPCRR when choosing different numbers of PCs.
- increase_prior_var.R: genomic prediction by BPCRR when manually varying the fixed prior variance.
- ridge_on_PCs.R: perform genomic prediction following a frequentist ridge regression on the PCs used for BPCRR, just as benchmark.

**NOTE on PC scaling**: In all the scripts, 'demo_sim.R', 'main.R', 'main_all_pis.R', 'main_accuracy_inspection.R' (UPDATE), when preparing the principal components, one can either scale PCs variances to 1 (standardize), by setting "scale = TRUE" in line 'XX <- scale(XX, center = TRUE, scale = FALSE)', or simply centering them keeping the line as it is, and later scaling them w.r.t. the variance of the first PC.

