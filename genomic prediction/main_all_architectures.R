

#### READ.ME ####

# Script for genomic prediction with BPCRR and BayesR
# Very similar to 'main.R', but handles a set of genetic architectures in one single run (it's just another loop).

rm(list=ls())
graphics.off()
setwd("/Users/zeliolina/Desktop/Università/tesi/code and data")

# Set the number of cores to use for parallel processing
options(mc.cores = 6)

#### Libraries
library(tictoc)
library(RSpectra)
library(Rfast)
library(INLA)
library(hibayes)


##### Functions
source("SVD.R")
source("BPCRR.R")
source("BPCRR_fix_effects.R")
source("BPCRR_PCeffect.R")
source("BPCRR_default.R")

append_row_and_save <- function(matrix_df, new_row, filename) {
  if (file.exists(filename)) {
    existing_data <- read.csv(filename, header = TRUE) # load
    # new_row <- as.data.frame(t(new_row))
    updated_matrix <- rbind(existing_data, new_row) # append
  }
  else {
    updated_matrix <- as.data.frame(t(new_row))
  }
  write.csv(updated_matrix, file = filename, row.names = FALSE) # write the updated matrix to the file
}

# Load a data frame from a file
load_matrix <- function(filename) {
  if (file.exists(filename)) {
    return(read.csv(filename))
  } else {
    return(NULL)
  }
}



#### Parameters settings ####

### number of SNPs to subsample
 n_subsample <- 10000
# n_subsample <- 50000
# n_subsample <- 100000

##### do.svd ####
# Set do.svd to 1 or 0 in order to do or not SVD of the current SNPs matrix. 
# [SVD is independent on the genetic architecture, thus we perform it once for each SNPs matrix (10k, 50k, 100k) and store the reduced matrices to retrieve them when needed]
do.svd = 0

##### Architectures ####
## E.g.
# a <- c(0.99, 0.95, 0.9, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05, 0.01)
# b <- c(0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99)

# specifying the sigmas
sigma1 <- 0
sigma2 <- 1e-4
sigma3 <- 1e-3
sigma4 <- 1e-2

### Vector of seeds
# we run a for loop where each time the seed is different so that we obtain distinct (but reproducible) values of accuracy
seeds <- c(788, 197823, 3456, 8765, 123)
len <- length(seeds)

## the following matrices are created only once and loaded each time this script is run.
## File names for saving the matrices for accuracy values
if (n_subsample == 10000) {
  BPCRR_file <- "BPCRR_default_accuracy_10K.csv"
  BayesR_file <- "BayesR_accuracy_10K.csv"
}

if (n_subsample == 50000) {
  BPCRR_file <- "BPCRR_accuracy_50K.csv"
  BayesR_file <- "BayesR_accuracy_50K.csv"
}

# 100k will be very very slow if you use up to 8 cores.
if (n_subsample == 100000) {
  # BPCRR_file <- "BPCRR_accuracy_100K.csv"
  BayesR_file <- "BayesR_accuracy_100K.csv"
  BPCRR_file <- "BPCRR_accuracy_100K_NEW.csv"
}


# Load existing matrices or create new empty ones
BPCRR <- load_matrix(BPCRR_file)
if (is.null(BPCRR)) {
  BPCRR <- data.frame(matrix(nrow = 0, ncol = len))
}
BayesR <- load_matrix(BayesR_file)
if (is.null(BayesR)) {
  BayesR <- data.frame(matrix(nrow = 0, ncol = len))
}



# ============================================================================ #


#### LOAD DATA ####

SNP <- readRDS("SNP_100k.rds")

# Dimensions of the SNP matrix as number of individuals and SNPs
dimensions <- dim(SNP)
Nanimals <- as.numeric(dimensions[1])
nSNPs <- n_subsample # as.numeric(dimensions[2])

# Store data as a matrix
SNP_unscaled <- as.matrix(SNP)


##### Subsample randomly
set.seed(36385)
if( n_subsample != 100000) {
  SNP_unscaled <- SNP_unscaled[, sample(ncol(SNP_unscaled), n_subsample)]
  SNP <- SNP[, sample(ncol(SNP), n_subsample)]
}


# ============================================================================ #


#### SVD ####

#### Load the already decomposed matrices you need
if (do.svd==0) {
  
  ## Load the saved SVD results you need
  if(n_subsample == 10000)
    load("svd_results_10K.RData")
  if(n_subsample == 50000)
    load("svd_results_50K.RData")
  if(n_subsample == 100000)
    load("svd_results_100K.RData")
  
  tmp.rspectra <- .GlobalEnv$svd_result
  
}


##### Perform SVD if not done yet
# Center column-wise, so that mean(col)=0, but do not scale to do the SVD
SNP_scaled <- scale(SNP,center=TRUE,scale=FALSE)

if (do.svd == 1) {
  source("SVD.R")
  
  # files to save the results
  if(n_subsample == 10000)
    filename_10K <- "svd_results_10K.RData"
  if(n_subsample == 50000 )
    filename_50K <- "svd_results_50K.RData"
  if(n_subsample==100000)
    filename_100K <- "svd_results_100K.RData"
  
  k <- 1200 # number of PCs to be estimated, default is the number of rows (3032)
  
  # SVD and save the results for each matrix ACCORDING TO THE NUMBER OF SNPS!
  if(n_subsample == 10000)
    tmp.rspectra <- perform_svd(SNP_scaled, k, filename_10K)
  if(n_subsample == 50000 )
    tmp.rspectra <- perform_svd(SNP_scaled, k, filename_50K)
  if(n_subsample==100000)
    tmp.rspectra <- perform_svd(SNP_scaled, k, filename_100K)
}

# Inspection the cumulative variance of the PCs
dd = prop.table(tmp.rspectra$d^2)
plot(cumsum(dd)) 


# ============================================================================ #

#### PCs ####

##### Ideal number of SNPs (formulas (5) to (7) in the manuscript)
h_k =0.33 * cumsum(dd)
N= 3032
k=seq(1:1200)

precision <- (N*h_k^2)/(N*h_k + k)
plot(k,precision)
which.max(precision) # 100k 779, 50k 773 #10k 763

# And the expected maximum precision
sqrt(precision[which.max(precision)]) # 100k 0.3931612, 50k 0.3933975

####
# This is (Odegards paper notation) T_q = X %*% V_q
# But now use the original matrix to multiply, not the one with core individuals!

# # to avoid compute it for 100k again
# XX  = SNP_scaled %*% tmp.rspectra$v
# file_name_XX <- "matrixXX.RData"
# save(XX, file = file_name_XX)
load("matrixXX.RData")

# Scaling all PCs with the standard deviation of the first PC
varPC1 <- var(XX[,1]) # 1
XX <- XX*((1/sqrt(varPC1)))

# And centering at the end of the preparation of the PCs (but not scaling the variances!)
XX <- scale(XX, center = TRUE, scale = FALSE)

# Preparing the PCs 
if( n_subsample == 10000)
  n_PCs <- 765   # PCs for 10k

if( n_subsample == 50000)
  n_PCs <- 775   #  773 PCs for 50k

if( n_subsample == 100000)
  n_PCs <- 780 # 779 PCs for 100k

XX_red <- XX[,1:n_PCs]

# Summing the total variance of the PCs to use as fixed prior
PCvar <- colVars(XX_red)

tot_PCvar <- sum(PCvar)

# As VA is set to 0.33, this goes into the fixed prior
varA <- 0.33
u.prior.var <- varA/tot_PCvar


# ============================================================================ #

#### Genomic prediction loop ####

for (k in 1:10) 
{
  # Determine the genetic architecture:
  pi1 <- a[k] # probability for an effect to be 0
  pi2 <- 0 # probability for an effect to be N(0,10^-4)
  pi3 <- 0 # probability for an effect to be N(0,10^-3)
  pi4 <- b[k] # probability for an effect to be N(0,10^-2)
  pi <- c(pi1, pi2, pi3, pi4)
 
  ## temporary vectors that will be appended to the respective matrices of results, at the end of the for loop
  v_accuracies_BPCRR <- c(0, 0, 0, 0, 0)  
  v_accuracies_BayesR <- c(0, 0, 0, 0, 0) 
  

  #### Inner loop: INLA and BayesR
  for (i in 1:len) {
    
    #### Simulate phenotypes using SNPs
    # randomly sampling integers from 0 to 3, with probabilities as specified by the pi's
    set.seed(seeds[i])
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
    
    #### Generate breeding values as LC of SNPs*effect sizes
    # Visualize effect sizes distribution
    hist(effect_sizes,nclass=50)
    # Generate breeding values as linear combinations of the SNPs times the effect sizes:
    aa <- SNP_unscaled %*% effect_sizes
    # Then scale the variance of the breeding values (Va) to 1
    aa <- aa/sd(aa)
    # Generate environmental effects where var = 2
    ee <- rnorm(Nanimals,0,sd=sqrt(2))
    # Summing breeding values and environmental effects into simulated phenotypes 
    yy <- aa + ee
    # Scaling so that var(yy) = 1, yielding a VA of 0.33
    yy <- scale(yy)
    
    #### Training-Test set
    # 5 folds such that each non-overlapping fold has approximately the same size
    # Randomly selecting the IDs in a vector, and thereby shuffling them
    # 20% circa of data as test
    set.seed(seeds[i])
    t.sample <- sample(1:Nanimals,round(Nanimals/5),replace=FALSE) #FALSE -> each index can only be selected once (unique subset)
    fold1 <- t.sample # sampled indices, which will be used as training set
    
    # Setting individual IDs
    IDC4 <- c(1:Nanimals)
    data <- as.data.frame(IDC4)
    # Adding the simulated phenotypes 
    data$yy <- yy
    data <- as.data.frame(data)
    
    # For this demo, we exclude the first fold from the fitting procedure
    # To this end, we are adding a new column in the data frame
    data$yy_CV <- data$yy
    
    # Replacing all phenotypes for the individuals in the first fold with NA. 
    # Note that INLA does then not use those rows to fit the model, which is exactly what we want. 
    # But we can then still get predictions for those entries
    fold <- as.vector(fold1)
    data$yy_CV[data$IDC4 %in% fold] <- NA
    
    # Implement BPCRR and BayesR, and save their results
    res_BPCRR <- BPCRR_fun(SNP_scaled, data, fold, Nanimals, yy_CV,u.prior.var)
    res_BayesR <- BayesR_fun(SNP_scaled, data, fold, Nanimals, yy_CV)
    v_accuracies_BPCRR[i] <- res_BPCRR[1]
    v_accuracies_BayesR[i] <- res_BayesR[1]

  }
  
  #### Update results ####
  
  # save to CSV files
  append_row_and_save(BPCRR,v_accuracies_BPCRR,BPCRR_file)
  append_row_and_save(BayesR,v_accuracies_BayesR,BayesR_file)
  
  }
  


#### Plot example ####
library(ggplot2)
library(dplyr)

# Define accuracy files and labels (choose notation and fix according to the architectures you used)
# labels <- c('(0.99 0 0 0.01)', '(0.95 0 0 0.05)', '(0.9 0 0 0.1)', '(0.8 0 0 0.2)', '(0.6 0 0 0.4)', '(0.4 0 0 0.6)', '(0.2 0 0 0.8)', '(0.1 0 0 0.9)', '(0.05 0 0 0.95)', '(0.01 0 0 0.99)')
labels <- c('0.99 | 0.01', '0.95 | 0.05', '0.9 | 0.1', '0.8 | 0.2', '0.6 | 0.4', '0.4 | 0.6', '0.2 | 0.8', '0.1 | 0.9', '0.05 | 0.95', '0.01 | 0.99')

accuracyBR <- read.table('BayesR_accuracy_100K.csv', header = TRUE, sep = ";")
accuracyBPCRR <- read.table('BPCRR_accuracy_scaleFALSE_100K.csv', header = TRUE, sep = ";")
accuracyBPCRR$Pi_s <- labels
accuracyBR$Pi_s <- labels

combined_df <- data.frame(
  Pi_s = rep(labels, each = 5),
  Accuracy = c(as.vector(t(accuracyBPCRR[,1:5])), 
               as.vector(t(accuracyBR[,1:5]))  ),
  Method = rep(c("BPCRR",
                 "BayesR"), 
               each = 5*length(labels)))


combined_df$Pi_s <- factor(combined_df$Pi_s, levels = labels)
combined_df$Method <- factor(combined_df$Method, 
                             levels = c("BPCRR", 
                                        "BayesR")) # Ensure order matches


ggplot(combined_df, aes(x = Pi_s, y = Accuracy, fill = Method)) +
  geom_boxplot() +
  labs(#title = "Genomic Prediction Comparison",
    x = "Genetic Architecture (π1 | π4)",
    y = "Accuracy") +
  scale_fill_manual(values = c("BPCRR" = "CornflowerBlue",
                               "BayesR" = "SkyBlue"  ),
                    labels = c("BPCRR",
                               "BayesR")
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    strip.background = element_blank(),  # Remove background for SNP facet labels
    strip.text = element_text(size = 10, face = "bold"),
    legend.title = element_blank(),  # Remove legend title
    legend.position = "bottom",
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(hjust = 0.5)
  ) +
  guides(fill = guide_legend(title = NULL))  # Simplify legend






