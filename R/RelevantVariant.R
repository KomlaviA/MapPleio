# res is the result from ComputePmS
RelevantVariant <-  function(res, treshold = 5*10^-8){
    IndList <- sapply(res, `[[`, 1) # table with subsets in columns and the indices of the traits in the corresponding subset
    PmList <- lapply(res, `[[`, 2) # list with Nbsubet elements, each element is a table with variants in rows and 3 columns (RSID, Pm, Pm_pvalue)
    tsub <- sapply(1:max(IndList), function(i) colSums(IndList == i)) # table of booleans, 1 if the trait in column is is in the subset in row, 0 otherwise
    PmPvalueMatrix <- sapply(PmList, `[[`, 3) # retrieve Pm_pvalues and convert into a matrix with variants in rows and subsets in column, at the intersection the p_values
    nbRelpvPerSNP <- outer(
    1 : nrow(PmPvalueMatrix),
    1 : ncol(tsub),
    Vectorize(function(i,j) length(which(PmPvalueMatrix[i, tsub[, j] == 1]<treshold))/length(which(tsub[,j] == 1)))
    )
    if (any(is.na(nbRelpvPerSNP)))
		warning("Some traits are not part of any subset.\n")
    return(nbRelpvPerSNP)
}



