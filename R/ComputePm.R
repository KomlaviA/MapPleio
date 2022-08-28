
ComputePm <- function(ZscoreWhitenedMatrix, RSids, LDCorrected = FALSE, POLYGENICITYCorrected = FALSE, NBSim = 25, GlobalTest = FALSE){
    # Compute the score
	Z0star2 <- 2
	Pm <- sqrt(rowSums(ZscoreWhitenedMatrix^2))

	# Load LD scores
	if(LDCorrected){
		cat("using LD scores from https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_v1.1_ldscores.tgz\n")
		load("Results/LDscores") # HOPS1 LD scores (from the server)
		warning("Matching LD scores on RSids, please make sure the SNPs are in build 37.\n")
		LDscores <- LDscores[match(RSids, LDscores$SNP), ]
		
		if(any(is.na(LDscores$baseL2)))
		warning(paste0(100 * sum(is.na(LDscores$baseL2)) / nrow(LDscores), "% of the SNPs will not be scored because of lacking LD scores.\n"))
		getLDCorrectedScore <- function(Pm, LDscores = LDscores){
			modPm <- lm(Pm ~ LDscores$baseL2)
			PmLD <- Pm - LDscores$baseL2 * modPm$coefficients[2]
			PmLD <- replace(PmLD, PmLD < 0, 0)
			return(PmLD)
		}
		lnormLD <- getLDCorrectedScore(Pm = Pm, LDscores = LDscores)
		Pm <- lnormLD # LD-corrected pleiotropy score (PmLD)
	}
    # Compute the p-values
	if(POLYGENICITYCorrected){
		getEmpirical <- function(ZscoreWhitenedMatrix, LD = LD){
			ZscoreWhitenedMatrixRandom <- apply(ZscoreWhitenedMatrix, 2, sample)
			lnormRandom <- cbind.data.frame(Pm = sqrt(rowSums(ZscoreWhitenedMatrixRandom^2)))
			if(LDCorrected)
				lnormRandom <- getLDCorrectedScore(Pm = lnormRandom$Pm, LDscores = LDscores)
			return(lnormRandom)
		}
        set.seed(123)
		EmpP <- replicate(NBSim, getEmpirical(ZscoreWhitenedMatrix, LD = LDCorrected))
		lnormRandom <- list(Pm = unlist(EmpP))
		Pm_percentiles <- rank(sort(lnormRandom$Pm) - 0.5)/length(lnormRandom$Pm)
		Pm_Pvalue <- 1 - Pm_percentiles[pmax(findInterval(Pm, sort(lnormRandom$Pm)),1)]

		Pm_Pvalue <- replace(Pm_Pvalue, Pm_Pvalue == 0, 0.99/(length(lnormRandom$Pm) * NBSim))

		if(GlobalTest){
			GT2 <- sapply(c("two.sided", "less", "greater"), function(alternative) wilcox.test(x = Pm^2, y = lnormRandom$Pm^2, exact = FALSE, alternative = alternative)$p.value)
		}

		# Empirical scores
		Pm <- sqrt(qchisq(p = Pm_Pvalue, df = ncol(ZscoreWhitenedMatrix), lower.tail = FALSE, log.p = FALSE))

	} else {
        # Theoretical Pvalues
		Pm_Pvalue <- pchisq(q = Pm^2, df = ncol(ZscoreWhitenedMatrix), lower.tail = FALSE, log.p = FALSE)

		if(GlobalTest){
			GT2 <- sapply(c("two.sided", "less", "greater"), function(alternative) wilcox.test(x = Pm^2, y = qchisq(p = seq(0, 1, length.out = length(Pm)), df = ncol(ZscoreWhitenedMatrix), lower.tail = FALSE, log.p = FALSE), exact = FALSE, alternative = alternative)$p.value)
		}
	}
    # Compile the results
	Scores <- cbind.data.frame(
		RSids = RSids,
		Pm = Pm / ncol(ZscoreWhitenedMatrix),
		Pm_Pvalue = Pm_Pvalue		
	)
	if(GlobalTest){
		return(list(cbind.data.frame(GlobalTest_Pm = GT2), Scores))
	} else {
		return(Scores)
	}
}
