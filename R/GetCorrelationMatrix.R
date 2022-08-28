library(FactoMineR)
library(bigstatsr)

# version 1 : produce some NaNs with th WLasso method
GetCorrelationMatrix <- function(ZscoreMatrix, method = c('WLasso','Empirical')){
	# Methods functions
	Empirical_Estimation <- function(ZscoreMatrix){
    	ZscoreFBM <- as_FBM(ZscoreMatrix)
		corBig <- big_cor(ZscoreFBM)
		corBig <- corBig[]
		return(corBig)
	}
	WLasso_Estimation <- function(ZscoreMatrix){
    	matrix_cor <- Empirical_Estimation(ZscoreMatrix) # compute the empirical correlation matrix (pxp)
    	CorApprox <- matrix(1, nrow(matrix_cor), ncol(matrix_cor))
    	groupsAuto <- HCPC(data.frame(matrix_cor), graph = FALSE, nb.clust = -1, consol = FALSE, method = "complete")$data.clust$clust
    					
    	for (i in which(table(groupsAuto) > 1)){ 
	   		matrixSubset <- matrix_cor[groupsAuto == i, groupsAuto == i]
	    	CorApprox[groupsAuto == i, groupsAuto == i] <- mean(matrixSubset[upper.tri(matrixSubset, diag = FALSE)])
    	}
    	# Inter-group correlation
    	AllComb <- combn(unique(groupsAuto), 2)
    	for (i in 1:ncol(AllComb)){
	    	matrixSubset <- matrix_cor[groupsAuto == AllComb[1, i], groupsAuto == AllComb[2, i]]
	    	# CorApprox[groupsAuto == AllComb[1, i], groupsAuto == AllComb[2, i]] <- mean(matrixSubset[upper.tri(matrixSubset, diag = FALSE)])
	    	# CorApprox[groupsAuto == AllComb[2, i], groupsAuto == AllComb[1, i]] <- mean(matrixSubset[upper.tri(matrixSubset, diag = FALSE)])
			CorApprox[groupsAuto == AllComb[1, i], groupsAuto == AllComb[2, i]] <- CorApprox[groupsAuto == AllComb[2, i], groupsAuto == AllComb[1, i]] <- mean(matrixSubset)
    	}
    	diag(CorApprox) <- 1
    	return(CorApprox)
	}
	method <- match.arg(method)
	# Correlation between the Zscores
	# If the method as argument is 'Wlasso', use the WLasso-like method to estimate the correlation matrix;
	# otherwise, use the empirical method.
	switch(method,  
	'WLasso'={
		ZscoreCorMatrix <- WLasso_Estimation(ZscoreMatrix)
    },
	'Empirical'={
		ZscoreCorMatrix <- Empirical_Estimation(ZscoreMatrix)
	})
	return(ZscoreCorMatrix)
}

