library(data.table)
library(bigsnpr)

# @param Variants: vector of variant IDs (chr:pos:A1:A2)
# @param PLINKdir: directory of the 1000Genomes project in bfile format
# @param plink: link to plink excecutable
# @param nbCores: number of cores to use for this function

# Returns a vector of independant variant IDs to keep after clumping

GetClumpedVariants <- function(Variants, PLINKdir = "/home/PleioMap/Data/1000Genomes/", plink = "/home/PleioMap/Scripts/plink", nbCores = 10, ...){

	CHR <- sort(unique(gsub("(^[0-9XYMT]+)(:.*)", "\\1", Variants)))
	
	SNPlist <- lapply(CHR, function(chr){
		cat(paste(chr, "\n"))
		bim <- fread(paste0(PLINKdir, "Chr", chr, ".bim"), data.table = FALSE)
		bim$variant56 <- paste(chr, bim$V4, bim$V5, bim$V6, sep = ":")
		bim$variant65 <- paste(chr, bim$V4, bim$V6, bim$V5, sep = ":")
		SNPs <- bim[bim$variant56 %in% Variants | bim$variant65 %in% Variants, 2]
		SNPs <- setdiff(SNPs, ".")
		write.table(SNPs, paste0(".tmpSNPs", chr) , row.names = FALSE, col.names = FALSE, quote = FALSE)
		cmd <- paste0(plink,
						" --bfile ", PLINKdir, "Chr", chr,
						" --extract .tmpSNPs",chr,
						" --threads ", nbCores, 
						" --make-bed --out .tmpChr", chr
		)
		system(cmd)
		obj.bed <- bed(paste0(".tmpChr", chr, ".bed"))
		KeptSNPs <- bed_clumping(obj.bed, ncores = 1, ...)
		KeptSNPs <- bim[bim[, 2] %in% SNPs[KeptSNPs], ]
		KeptSNPs <- Variants[Variants %in% KeptSNPs$variant56 | Variants %in% KeptSNPs$variant65]
	})
	system("rm .tmp*")
	SNPlist <- unlist(SNPlist)
	return(SNPlist)
}
