# In this script, we demonstrate how estimation of VA and genomic prediction
# was performed with BPCRR using fixed priors (fomula (4) in the manuscript) and BayesR
# For the purpose of this demonstration, only 1000 SNPs are used. Potentially this is too little.

# Set to local directory
setwd("~/Subversion/PhDJanne/demo_sim")

### Necessary packages:
library(tictoc)
library(RSpectra)
library(Rfast)
library(INLA)
library(hibayes)

### Setting some variables
# Should svd be run? This parameter is 1 if you run the script for the first time, but may be set to 0 afterwards, so you avoid doing the svd calculations repeatedly.
do.svd=1
# Set a seed (you may try various ones)
seed <- 788

# Loading the SNPs (3032 individuals with 1000 SNPs each - maybe these are too few and can be expanded later)
#load("~/Subversion/PhDJanne/demo_sim/SNPs_demo.Rdata")
#d.SNP <- SNP_SNP

#d.SNP <- readRDS("~/Subversion/PhDJanne/demo_sim/SNP_10k.rds")
d.SNP <- readRDS("~/Subversion/PhDJanne/demo_sim/SNP_50k.rds")
#d.SNP <- readRDS("~/Subversion/PhDJanne/demo_sim/SNP_100k.rds")



# Dimensions of the SNP matrix as number of individuals and SNPs
(dimensions <- dim(d.SNP))
Nanimals <- as.numeric(dimensions[1])
nSNPs <- as.numeric(dimensions[2])

# Store the data as a matrix
d.SNP_unscaled <- as.matrix(d.SNP)


#############   Using SNPs to generate simulated phenotypes
# Determine the genetic architecture:
pi1 <- 0.95 # probability for an effect to be 0
pi2 <- 0.04   # probability for an effect to be N(0,10^-4)
pi3 <- 0.008   # probability for an effect to be N(0,10^-3)
pi4 <- 0.002   # probability for an effect to be N(0,10^-2)

(pi1+pi2+pi3+pi4)

# specifying the sigmas
sigma1 <- 0
sigma2 <- 1e-4
sigma3 <- 1e-3
sigma4 <- 1e-2

# randomly sampling integers fro 0 to 3, whit probabilities as specified by the pi's.
set.seed(seed)
dist_choice <- sample(c(0:3), nSNPs, replace = TRUE, prob = c(pi1,pi2,pi3,pi4))

effect_sizes <- c(NA) 

# setting the effect sizes for each SNP
for (ii in 1:nSNPs){
  if (dist_choice[ii]==0) {
    effect_sizes[ii] <- rnorm(1,0,sqrt(sigma1))
  }
  if (dist_choice[ii]==1) {
    effect_sizes[ii] <- rnorm(1,0,sqrt(sigma2))
  }
  if (dist_choice[ii]==2) {
    effect_sizes[ii] <- rnorm(1,0,sqrt(sigma3))
  }
  if (dist_choice[ii]==3) {
    effect_sizes[ii] <- rnorm(1,0,sqrt(sigma4))
  }
}


# Visually checking the distribution of effect sizes.
hist(effect_sizes,nclass=50)

# Generate breeding values as linear combinations of the SNPs times the effect sizes:
aa <- d.SNP_unscaled %*% effect_sizes
# Then scale the variance of the breeding values (Va) to 1
aa <- aa/sd(aa)

# Generate environmental effects where var = 2
ee <- rnorm(Nanimals,0,sd=sqrt(2))

# Summing breeding values and environmental effects into simulated phenotypes 
yy <- aa + ee

# Scaling so that var(yy) = 1, yielding a VA of 0.33
yy <- scale(yy)

#######
# Here we are generating 5 folds such that each non-overlapping fold has approximately the same size.
# Randomly selecting the IDs in a vector, and thereby shuffling them
set.seed(seed)
t.sample <- sample(1:Nanimals,round(Nanimals/5),replace=FALSE)

fold1 <- t.sample
  
# 
# # And specifying the 10 folds
# fold1 <- Shuffle[1:606]
# fold2 <- Shuffle[607:1212]
# fold3 <- Shuffle[1213:1819]
# fold4 <- Shuffle[1820:2425]
# fold5 <- Shuffle[2426:3032]



#######. Performing the SVD
# Note that the SVD is totally independent on the genetic architecture, thus the expected max precision is also fixed for a given 

# Center column-wise, so that mean(col)=0, but do not scale to do the SVD
d.SNP_scaled <- scale(d.SNP,center=TRUE,scale=FALSE)


if (do.svd==1){
  k <- 1200 # number of PCs to be estimated, default is # rows
  
  # Use the svds() function from the RSpectra package
  # Max iterations for our experiments was set to 2*10^6, instead of the default 1*10^6
  tmp.rspectra <- svds(d.SNP_scaled, k, nu = k, nv = k, opts = list(maxitr = 2000000))
  
  # Inspection the cumulative variance op the PCs
  dd = prop.table(tmp.rspectra$d^2)
  plot(cumsum(dd))
}


###### Look at the formula (formulas (5) to (7) in the manuscript) to find ideal number of SNPs
h_k =0.33 * cumsum(dd)
N= 3032
k=seq(1:1200)

precision <- (N*h_k^2)/(N*h_k + k)

plot(k,precision)

which.max(precision)
# Shows that roughly 770 PCs are optimal. So let's go for it.

# And the expected maximum precision is the given as
sqrt(precision[which.max(precision)])


######
# This is (in Odegards paper notation) T_q = X %*% V_q
XX  = d.SNP_scaled %*% tmp.rspectra$v

# Scaling all PCs with the standard deviation of the first PC
varPC1 <- var(XX[,1])
XX <- XX*((1/sqrt(varPC1)))

# And centering (but not scaling the variances!) at the end of the preparation of the PCs
XX <- scale(XX, center = TRUE, scale = FALSE)

# Preparing the PCs 
n_PCs <- 800 # Deciding to use 550 according to the above formula
XX_red <- XX[,1:n_PCs]


# Summing the total variance of the PCs to use as fixed prior (formula (4))

PCvar <- colVars(XX_red)
tot_PCvar <- sum(PCvar)

# As VA is set to 0.33, this goes into the fixed prior
varA <- 0.33
u.prior.var <- varA/tot_PCvar


# Setting up the folds for cross validation, (stone age style!), to be able to inspect all results on a individual level over all traits, both simulated and real traits

# Setting individual IDs
IDC4 <- c(1:Nanimals)
data <- as.data.frame(IDC4)
# Adding the simulated phenotypes 
data$yy <- yy
data <- as.data.frame(data)


# For this demo, we exclude the first fold from the fitting procedure
# To this end, we are adding a new column in the data frame
data$yy_CV <- data$yy

# Replacing all phenotypes for the individuals in the first fold with NA. Note that INLA does then not use those rows to fit the model, which is exactly what we want. But we can then still get predictions for those entries
fold <- as.vector(fold1)
data$yy_CV[data$IDC4 %in% fold] <- NA


############################################
### INLA setup fixed prior
############################################

# In inla, we use the latent z-model to model the marker effects with ridge shrinkage. The formula is given as
formula.sim  = yy_CV ~ f(IDC4, model = "z", Z = XX_red,
                     hyper=list(
                       prec=list(initial=log(1/u.prior.var),fixed=TRUE) # fixed=FALSE or fixed=TRUE determines whether the variance is fixed or not. 
                     ))

# Actually calling inla() - this call takes a bit of time:
model.sim = inla(formula=formula.sim, family="gaussian",
                  data=data,
                  control.family=list(hyper = list(theta = list(initial=log(3/2),  
                                                                prior="pc.prec",
                                                                param=c(1,0.1)))), 
                  control.compute=list(dic=F, config = FALSE, 
                                       return.marginals=FALSE
                                       ), # To be able to resample from the inla object, we need config=TRUE. However, config=FALSE makes the computation faster, so only set it to TRUE if you plan to resample (i.e., if you want to estimate VA)
                  num.threads=10 # using only 10 out of the 12 cores I have
)

#### Summary for inla objects
summary(model.sim)


### Evaluate prediction accuracy ###
# Finding the breeding values and correlation between predicted breeding value and the simulated phenotype
breedingv <- model.sim$summary.random$IDC4$mean[1:Nanimals]

IDC4 <- c(1:Nanimals)

# Data frame with IDs and breeding values
bvAndIds <- as.data.frame(IDC4)
bvAndIds$breedingv <- breedingv

# Data frame with both old and predicted phenotypes as well as breeding values
result_allIDs <- merge(x=data,y=bvAndIds,by="IDC4",all.x=TRUE)

# Data frame with both old and predicted phenotypes as well as breeding values for the IDs in the respective fold, and plotting the old phenotype against the breeding values
result_foldIDs <- result_allIDs[result_allIDs$IDC4 %in% fold,]
plot(result_foldIDs$yy, result_foldIDs$breedingv)

# And finally the correlation between the predicted breeding value and the actual phenotypes that were spared out in the respective fold:
(cor_phenobreed_sim <- cor(result_foldIDs$yy, result_foldIDs$breedingv,use="complete.obs"))

# Correlation between true and predicted breeding values
cor_phenobreed_sim

##### Resampling from the posteriors - if interest centers around variance components.
# nsamples <- 500
# sample.sim <- inla.posterior.sample(n=nsamples,model.sim)
# 
# # Extract samples of the breeding values; in each entry of the list the breeding values of all animals in one sample are stored
# samples.a <-  lapply(1:nsamples, function(ii) {
#   sample.sim[[ii]]$latent[substring(rownames(sample.sim[[1]]$latent),1,4)=="IDC4"][1:Nanimals]     
# })
# 
# VAs.sim <- unlist(lapply(samples.a, var))
# mean(VAs.sim) # Estimated VA
# hist(VAs.sim)



###################################
####  BayesR setup and computations
###################################

M <- as.matrix(d.SNP_scaled)  # SNP matrix

M.id <- c(1:3032) # Vector of IDs
 
fit <- ibrm(yy_CV ~ 1, 
            data = data,
            M = M,
            M.id = M.id,
            method = "BayesR", 
            niter=10000, nburn=5000,
            thin = 10)

# Finding the breeding values
breedingv <- fit$g
breedingv <- breedingv$gebv

IDC4 <- c(1:Nanimals)

# Data frame with IDs and breeding values
bvAndIds <- as.data.frame(IDC4)
bvAndIds$breedingv <- breedingv

# Data frame with both old and predicted phenotypes as well as breeding values
result_allIDs <- merge(x=data,y=bvAndIds,by="IDC4",all.x=TRUE)

# Data frame with both old and predicted phenotypes as well as breeding values for the IDs in the respective fold, and plotting the old phenotype against the breeding values
result_foldIDs <- result_allIDs[result_allIDs$IDC4 %in% fold,]
plot(result_foldIDs$yy, result_foldIDs$breedingv)

# And finally the correlation between the presicted breedingvalye and the 
# old phenotype
(cor_phenobreed_bayesR <- cor(result_foldIDs$yy, result_foldIDs$breedingv,use="complete.obs"))

cor_phenobreed_bayesR
