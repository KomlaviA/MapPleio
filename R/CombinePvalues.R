CombinePvalues <-  function(PmScores, NbTraitSubset) {
    PmList <- lapply(PmScores, `[[`, 2)
    PmMatrix <- sapply(PmList, `[[`, 2) # convert into a matrix with variants in rows and subsets in column, at the intersection the score Pm
    MeanPm <- sqrt(rowSums(PmMatrix^2))/ncol(PmMatrix) # average score
    PmPval <- sapply(PmList, `[[`, 3)
    # Pvcomb <- sapply(1:nrow(PmPval), function(i) pchisq(-2 * sum(log(i)), df=2*length(PmList), log.p = TRUE, lower.tail=FALSE))
    Statistic <- apply(PmPval, MARGIN = 1, FUN = function(pvals) max(pvals))
    # S <- replace(S, S == 0, min(setdiff(S, 0)))
    Statistic <- as.matrix(Statistic)
      
    # p-value
    P_value <- Statistic
    P_value <- replace(P_value, P_value == 0, min(setdiff(P_value, 0)))
        
    # Pm score
    Pm <- qchisq(p = P_value, df = NbTraitSubset, lower.tail = FALSE, log.p = FALSE)
    
    comp <- cbind(MeanPm, P_value, Pm)
    colnames(comp) <- c("Average_Pm","Combined_Pvalue ", "Combined_Pm")
    # Pmcomb <- sqrt(qchisq(p = Pvcomb, df = NbTraitSubset, lower.tail = FALSE, log.p = FALSE))/NbTraitSubset
    return(comp)
}
      
