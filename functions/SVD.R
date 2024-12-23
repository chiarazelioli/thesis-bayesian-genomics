perform_svd <- function(matrix, k, filename) {
  svd_result <- svds(matrix, k, nu = k, nv = k, opts = list(maxitr = 2000000))
  save(svd_result, file = filename)   # save the SVD 
  return(svd_result)
}
