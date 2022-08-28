#' Title
#'
#' @param res Results from the ComputePmS function
#'
#' @return An array with traits in rows and three columns: Mean with traits, Mean without traits, Pvalue
#' @export
#'
#' @examples
#' Trimp <- TraitInfluence(Pmsamp)
TraitInfluence <- function(res){
  IndList <- sapply(res, `[[`, 1) # table with subsets in columns and the indices of the traits in the corresponding subset
  PmList <- lapply(res, `[[`, 2) # list with Nbsubet elements, each element is a table with variants in rows and 3 columns (RSID, Pm, Pm_pvalue)
  tsub <- sapply(1:max(IndList), function(i) colSums(IndList == i)) # table of booleans, 1 if the trait in column is is in the subset in row, 0 otherwise
  PmMatrix <- sapply(PmList, `[[`, 2) # convert into a matrix with variants in rows and subsets in column, at the intersection the score Pm
  Means <- function(i, PmMatrix, tsub){
    Indices <- tsub[,i]
    if(sum(Indices) == 1)
      MeanW <- PmMatrix[, Indices == 1]
    else
      MeanW <- apply(PmMatrix[, Indices == 1], MARGIN = 1, FUN = mean) # mean per variant of the scores from the subsets with the trait # MARGIN = 1 (rows; 2 is for columns); FUN = mean (mean of Pm from subsets with the trait)
    MeanO <- apply(PmMatrix[, Indices == 0], MARGIN = 1, FUN = mean) # mean per variant of the scores from the subsets without the trait
    if (any(is.na(MeanW)) | any(is.na(MeanO)))
      p_value <- NA
    else
      p_value <- t.test(MeanW, MeanO, paired = FALSE)$p.value
    meanres <- c(mean(MeanW), mean(MeanO), p_value)
    return(meanres)
  }
  TraitPm <- t(sapply(1:ncol(tsub), Means, PmMatrix = PmMatrix, tsub = tsub))
  colnames(TraitPm) <- c("Mean with trait", "Mean without trait", "Pvalue")
  if (any(is.na(TraitPm)))
    warning(paste("Some traits have not been evaluated, the number of subsets of", length(IndList), "might not be high enough.\n"))
  return(TraitPm)
}

# The results from ComputePmS contains subsets with SNPs and traits and their Pm scores and Pm_pvalue
# The function determinate the importance of traits by computing the mean of Pm scores in the two cases : with and without trait.
# retrieve the subsets in which each feature is found
