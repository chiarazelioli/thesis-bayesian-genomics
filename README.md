# Master Thesis - work in progress
Benchmarking existing state-of-the-art Bayesian methods (BPCRR and BayesR) for genomic prediction, investigating their performance across different genetic architectures.

TODO: Explain or link references for BPCRR and BayesR.
      Explain step by step the scripts when you organise the all repository and insert the scripts.

The dataset, "SNP_100k.rds", contains 100'000 SNPs for 3032 genotyped individuals. Contact stefamu@ntnu.no to receive the link to download it.

In 'demo_sim.R' we demonstrate how estimation of VA and genomic prediction are performed via BPCRR using fixed priors and BayesR. A user has to set (explain...)

Refer to 'main.R' to perform genetic prediction running several (5, but one can potentially increase it by setting more seeds) simulations for the same genetic architecture π that the user has to set, in order to obtain several instances of genomic prediction accuracy and visually compare through a box-plot the performance of the two methods BPCRR and BayesR. This script optimises several procedures, and makes use of two external functions to apply each method: 'BPCRR.R' and 'BayesR.R'. 'BPCRR_default.R' implements BPCRR using the default, uninformative, Gamma prior. 'SVD.R' is the function that is called to perform singular value decomposition of the SNP matrix whenever the binary parameter 'do.svd' is set to 1.

To run all genetic architectures at once, refer to 'main_all_pis.R', which contains an outer loop that iterates over several genetic architectures πs, simulating from an oligogenic to a highly polygenic trait. 

'Ridge.S'

'main_accuracy_inspection' FIX THE NAME AND SPLIT THE SCRIPT

TODO: another BPCRR function to account for fixed effects

NOTE on PC scaling: In all the scripts, 'demo_sim.R', 'main.R', 'main_all_pis.R', 'main_accuracy_inspection.R' (UPDATE), when preparing the principal components, one can either scale PCs variances to 1 (standardize), by setting "scale = TRUE" in line 'XX <- scale(XX, center = TRUE, scale = FALSE)', or simply centering them keeping the line as it is, and later scaling them w.r.t. the variance of the first PC.

