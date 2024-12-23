#### BayesR ####
BayesR_fun <- function(d.SNP_scaled, data, fold, Nanimals, yy_CV) {
  
  tic()
  
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
  
  # And finally the correlation between the predicted breeding value and the old phenotype
  (cor_phenobreed_bayesR <- cor(result_foldIDs$yy, result_foldIDs$breedingv,use="complete.obs"))
  
  BayesR_time <- toc()
  
  res <- c(cor_phenobreed_bayesR, BayesR_time$toc - BayesR_time$tic)
  return(res)
}
