#### INLA - fixed prior ####

BPCRR_fun <- function(SNP_scaled, data, fold, Nanimals, yy_CV, u.prior.var) {
  

tic()

# Use the latent z-model to model the marker effects with ridge shrinkage
formula.sim  = yy_CV ~ f(IDC4, model = "z", Z = XX_red,
                         hyper=list(
                           prec=list(initial=log(1/u.prior.var),fixed=TRUE) # fixed=FALSE or fixed=TRUE determines whether the variance is fixed or not.
                         ))


# Call inla()
model.sim = inla(formula=formula.sim, family="gaussian",
                 data=data, # also this you should add the 5 cols at the beginning of the dataframe
                 control.family=list(hyper = list(theta = list(initial=log(3/2),  
                                                               prior="pc.prec",
                                                               param=c(1,0.1)))), 
                 control.compute=list(dic=F, config = FALSE, 
                                      return.marginals=FALSE
                 ), # To be able to resample afterwards, from the INLA object, we need config=TRUE. However, config=FALSE makes the computation faster, so only set it to TRUE if you plan to resample (i.e., if you want to estimate VA)
                 
                 # num.threads=10 # number of chores
)


summary(model.sim)


## Evaluate prediction accuracy
# by finding the breeding values and correlation between predicted breeding value and the simulated phenotype
breedingv <- model.sim$summary.random$IDC4$mean[1:Nanimals]

IDC4 <- c(1:Nanimals)

# Data frame with IDs and corresponding breeding values
bvAndIds <- as.data.frame(IDC4)
bvAndIds$breedingv <- breedingv

# Data frame with both old and predicted phenotypes as well as breeding values
# merge the original dataframe 'data' with a new dataframe 'bvAndIds'  
result_allIDs <- merge(x=data,y=bvAndIds,by="IDC4",all.x=TRUE)

# Data frame with both old and predicted phenotypes as well as breeding values for the IDs in the respective fold, and plotting the old phenotype against the breeding values
# filters the merged data frame result_allIDs to include only the rows where the IDC4 is in the current fold (fold). 
# This creates result_foldIDs, which contains only the subset of data corresponding to the current fold being processed
result_foldIDs <- result_allIDs[result_allIDs$IDC4 %in% fold,]
# plot(result_foldIDs$yy, result_foldIDs$breedingv)

# And finally the correlation between the predicted breeding value and the actual phenotypes that were spared out in the respective fold:
(cor_phenobreed_sim <- cor(result_foldIDs$yy, result_foldIDs$breedingv,use="complete.obs"))


##### Resampling from the posteriors #####
## if interested, centers around variance components. Go up and set "config=TRUE" when fitting the model with INLA

# nsamples <- 500
# sample.sim <- inla.posterior.sample(n=nsamples,model.sim)

# # Extract samples of the breeding values; in each entry of the list the breeding values of all animals in one sample are stored
# samples.a <-  lapply(1:nsamples, function(ii) {
#   sample.sim[[ii]]$latent[substring(rownames(sample.sim[[1]]$latent),1,4)=="IDC4"][1:Nanimals]  
# })

# VAs.sim <- unlist(lapply(samples.a, var))

## Estimated VA
# mean(VAs.sim) 
# hist(VAs.sim)

INLA_time <- toc()

# if in need, return VA too
res <- c(cor_phenobreed_sim, INLA_time$toc - INLA_time$tic)
return(res)

}




