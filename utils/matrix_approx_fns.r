# A collection of functions that perform or help perform matrix approximations
#   (saved here to help the structure of code)
#   (often repeated functions/purposes from different sims/analysis- moved some fns here)
# These functions should have the general inputs:
# 
# Overall Main Function:
#   - approx_matrix- will call one of the following based on the specified method
# 
# These methods include:
#   - SVD
#   - Sparse SVD
#   - softImpute
#   - spectralbiclust
#   - zeros
#   - average
#
#


require(softImpute)



#' Perform Matrix Completion for leave one out X
#' Previously in matrixPrior_utils.R. 
#' Matrix approx functions that implement softImpute approximations require this function 
#' (e.g. matrix_shrinkage.r's approx_matrix(...))
#' '
#' @param X (matrix)
#' @param type (character) 'svd' or 'als' for the type of softImpute method used
create_matrix_completion_softImpute <- function(X, rank_max, type = 'svd') {
  # debug
  # X=matrix(rnorm(30),6,5)
  # X[sample(1:30,10,replace=FALSE)]=NA
  
  
  
  
  
  # Part 1: Set up hyperparameters/initial start to reuse ?? 
  # choosing lambda parameter ??? idk, just do (1/2) * lam0 for now??
  X_dgCMatrix=as(X,"Incomplete")
  lam0=lambda0(X_dgCMatrix)
  # 
  # lamseq=exp(seq(from=log(lam0),to=log(1),length=10))
  # 
  # fits=as.list(lamseq)
  # ranks=as.integer(lamseq)
  # rank.max=2
  # warm=NULL
  # nuclear_norms = as.list(lamseq)
  # matrix_completed = as.list(lamseq)
  # for( i in seq(along=lamseq)){
  #   fiti=softImpute::softImpute(X_dgCMatrix,lambda=lamseq[i],rank=rank.max,warm=warm)
  #   ranks[i]=sum(round(fiti$d,4)>0)
  #   rank.max=min(ranks[i]+2,4)
  #   warm=fiti
  #   fits[[i]]=fiti
  #   # cat(i,"lambda=",lamseq[i],"rank.max",rank.max,"rank",ranks[i],"\n")
  #   matrix_completed[i] = softImpute::complete(X,fiti,unscale=TRUE)
  #   nuclear_norms[i] = sum(fiti$d)
  # }
  # 
  # # complete matrix
  # X_completed = softImpute::complete(X,fiti,unscale=TRUE)
  # X_completed = matrix(NA, nrow=nrow(X), ncol=ncol(X)) 
  # for(i in 1:nrow(X)) {
  #   for(j in 1:ncol(X)) {
  #     
  #   }
  # }
  # 
  
  
  
  #' @param i (integer) row index
  #' @param j (integer) col index
  matrix_completion_softImpute <- function(i, j) {
    X_ij = X
    X_ij[i,j] = NA
    # use the same fit for i,j element complete (used for sample splitting version)
    fit = softImpute::softImpute(X_ij, lambda = (lam0 + 1)/2, rank=rank_max, warm=NULL, type = type)
    
    
    X_completed = softImpute::complete(X_ij, fit, unscale=TRUE)
    
    return(X_completed[i,j])
  }
  
  return(matrix_completion_softImpute)
  
  
}




