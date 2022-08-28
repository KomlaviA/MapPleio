#' Compute the pleiotropy magnitude score (Pm) of SNPs in a subset of traits
#'
#' @param ZscoreWhitenedMatrix The whitened Zscore matrix from GetWhitenedZscores function.
#' @param NbTraitSubset The number of traits in a subset
#' @param NbSubset The number subset
#' @param ... others parameters from the ComputePm function
#'
#' @return
#' @export
#'
#' @examples
#' Pmperssample <- ComputePmS(ZscoreWhitenedMatrix = ZscoreMatrixWhitened, RSids = row.names(ZscoreMatrixWhitened), LDCorrected = FALSE, POLYGENICITYCorrected = FALSE, NBSim = 25, GlobalTest = FALSE, NbSubset = 10)

ComputePmS <- function(ZscoreWhitenedMatrix, NbTraitSubset = round(ncol(ZscoreWhitenedMatrix) / 3), NbSubset = 10, ...){
  GetSubsets <- function(i){
    TraitIndices <- sample(1:ncol(ZscoreWhitenedMatrix), size = NbTraitSubset, replace = FALSE) # get the indices of traits
    ZscoreWhitenedMatrix_subsets <- ZscoreWhitenedMatrix[, TraitIndices] # get the whitened Z score matrix
    Pm <- ComputePm(ZscoreWhitenedMatrix_subsets, ...) # ... allows to take into account all the parameters from function Compute_Pm()
    PmList <- list(TraitIndices, Pm)
    return(PmList)
  }
  res <- lapply(1:NbSubset, FUN = GetSubsets)
  return(res)
}


# GetSubsets
# This function takes as in input a given subset of traits in columns and SNPs in rows.
# It computes the Pm score for each SNPs in each subset.
