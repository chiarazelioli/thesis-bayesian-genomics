


# ############################################
# ### INLA setup default prior
# ############################################

BPCRR_default_fun <- function(SNP_scaled, data, fold, Nanimals, yy_CV, u.prior.var) {
  
  formula  = yy_CV ~ f(IDC4, model = "z", Z = XX_red,
                       hyper=list(
                         prec=list(initial=log(1/u.prior.var),fixed=FALSE)
                       ))
  
  model.sim.default = inla(formula=formula, family="gaussian",
                           data=data,
                           control.family=list(hyper = list(theta = list(initial=log(1),  #1
                                                                         prior="pc.prec",
                                                                         param=c(1,0.1)))), 
                           control.compute=list(dic=F, config = TRUE, return.marginals=FALSE) #to be able to resample afterwards
  )
  
  
  
  # #### Summary
  # summary(model.sim.default)
  # 
  # #### Resampling from the posteriors of Va
  # nsamples <- 500
  # sample.sim.default <- inla.posterior.sample(n=nsamples,model.sim.default)
  # 
  # # To extract samples of the breeding values; in each entry of the list the breeding values of all animals in one sample are stored
  # samples.a <-  lapply(1:nsamples, function(ii) {
  #   sample.sim.default[[ii]]$latent[substring(rownames(sample.sim.default[[1]]$latent),1,4)=="IDC4"][1:Nanimals]
  # })
  # 
  # VAs.sim.default <- unlist(lapply(samples.a, var))
  # mean(VAs.sim.default) # Estimated VA
  # hist(VAs.sim.default)
  
  
  ### Evaluate prediction accuracy ###
  # Finding the breeding values and correlarion between predicted breedingvalue
  # and the simulated phenotype
  breedingv <- model.sim.default$summary.random$IDC4$mean[1:Nanimals]
  
  IDC4 <- c(1:Nanimals)
  
  # Dataframe with IDs and breedingvalues
  bvAndIds <- as.data.frame(IDC4)
  bvAndIds$breedingv <- breedingv
  
  # Dataframe with both old and predicted phenotypes as well as breedingvalues
  result_allIDs <- merge(x=data,y=bvAndIds,by="IDC4",all.x=TRUE)
  
  # Dataframe with both old and predicted phenotypes as well as breedingvalues
  # for the IDs in the respective fold, and plotting the old phenotype
  # against the breedingvalues
  result_foldIDs <- result_allIDs[result_allIDs$IDC4 %in% fold,]
  # plot(result_foldIDs$yy, result_foldIDs$breedingv)
  
  # And finally the correlation between the predicted breeding value and the actual phenotypes that were spared out in the respective fold:
  (cor_phenobreed_sim_default <- cor(result_foldIDs$yy, result_foldIDs$breedingv,use="complete.obs"))
  
  return(cor_phenobreed_sim_default)

}







