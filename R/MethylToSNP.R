#' A novel tool for removing SNPs in Illumina DNA methylation array data
#'
#' @docType package
#' @name MethylToSNP
NULL


# require(Ckmeans.1d.dp)

#' Identify sites that may have underlying SNPs in methylation array data
#' 
#' @param data An GenomicRatioSet, GenomicMethylSet, MethylSet, or RatioSet object (see minfi package)
#' @param gap.ratio The ratio of two gaps should be above the threshold.
#' @param gap.sum.ratio The ratio of the sum of two gaps relative to the total range of values should be above the threshold.
#' @param outlier.sd Do not consider outliers that are more than the specified number of standard deviations from the cluster center.
#' @param SNP Optional SNP annotation
#' @param verbose Show additional information. Useful for debugging.
#' @return Detected probes with 3-tier SNP-like methylation pattern along with their reliability scores and SNP annotation
#' @examples
#' MethylToSNP(data)
#' @export
MethylToSNP <- function(data, gap.ratio = 0.75, gap.sum.ratio = 0.5, verbose=FALSE, outlier.sd = 3.0, SNP)
{
    if ( (!length(dim(data)) == 2) || (dim(data)[2] < 2)){
        stop("[MethylToSNP] There have to be at least 2 samples")
    }

    if (dim(data)[2] < 50){
        message("[MethylToSNP] Warning, SNP detection may be unreliable in datasets with less than 50 samples")
    }
        
	if (gap.ratio >= 1 || gap.ratio <= 0){
		stop("[MethylToSNP] You have entered an unacceptable gap.ratio value. This must be a number between zero and one.")
	}
	
    if (gap.sum.ratio <= 0 || gap.sum.ratio >= 1){
		stop("[MethylToSNP] You have entered an unacceptable gap.sum.ratio value. This must be a number between zero and one.")
	}
    
    if (is(data, "GenomicRatioSet") || is(data, "GenomicMethylSet") || is(data, "MethylSet") || is(data, "RatioSet")) {
        if (verbose) 
            message("[MethylToSNP] Calculating beta matrix. minfi is required")
        .checkMinfi()
        data <- getBeta(data)
    }

    if (missing(SNP)){
    	message("[MethylToSNP] Optionally, specify SNPs in a data frame with row names corresponding to cg probes (such as SNPs.147CommonSingle in minfiData or minfiDataEPIC package)")
    }

    ##########
    # STEP 1.
    # Identify probes that could be potential SNPs in a given set of samples
    
	potentialSNPs <- NULL
	probes <- rownames(data)

	# for each probe:
	for (i in 1: dim(data)[1]){

		probe <- probes[i]
		x <- as.numeric(na.omit(data[i,]))
		if (length(x) < 3) {
			next
		}
		
		span <- diff(range(x, na.rm = TRUE))
		if (span >= 0.5) {
			
			# inverse density of beta values
			weights <- 1.0 / approxfun(density(x))(x)

			# optimal 1D clustering with dynamic programming into 3 clusters
			# no randomization involved
			kmeans <- Ckmeans.1d.dp(x=t(x), y=weights, k=3)
			clusters <- kmeans$cluster
			centers <- kmeans$centers

			top <- which(centers == max(centers))
			bottom <- which(centers == min(centers))
			middle <- c(1:3)[-c(top, bottom)]


			# disregard outliers (assign them to non-existing cluster #4)
			if (outlier.sd != FALSE) {
				top_outliers <- which(abs(scale(x[clusters == top])) > outlier.sd)
				clusters[top_outliers] <- 4

				middle_outliers <- which(abs(scale(x[clusters == middle])) > outlier.sd)
				clusters[middle_outliers] <- 4

				bottom_outliers <- which(abs(scale(x[clusters == bottom])) > outlier.sd)
				clusters[bottom_outliers] <- 4

				# if we removed outliers and one of the clusters is empty -- should not happen
				if (length(x[clusters == top]) * length(x[clusters == middle]) * length(x[clusters == bottom]) == 0) {
					next
				}
			}

			# two gaps:
			# min of the top cluster minus max of the middle cluster
			# and min of the middle cluster minus max of the bottom cluster			
			gaps.sorted <- sort(c(
				min(x[clusters == top]) - max(x[clusters == middle]),
				min(x[clusters == middle]) - max(x[clusters == bottom])),
				decreasing = TRUE)

			gaps.largest <- gaps.sorted[1]
			gaps.smallest <- gaps.sorted[2]
					
			###
		    # Apply gap thresholds to decide whether to add a cg probe to the list of potential SNPs
		    #
			if ((gaps.smallest >= gap.ratio * gaps.largest) && (sum(gaps.sorted) >= gap.sum.ratio * span)) {
			    if (verbose) {
			        message(probe)
			    }
			    potentialSNPs <- append(potentialSNPs, probe)
			}
		} else {
			# data span (max - min) is too narrow
		}
		
		###
		# Progress indicator
		#
		if (verbose){
			if((i %% 1000) == 0){
				message('[MethylToSNP] Processed: ', i, " Identified potential SNPs: ", length(potentialSNPs))
			}
		}
	}
	###
	# Summary
	#
	if(verbose){
		if (length(potentialSNPs) > 0 ){
		    message("[MethylToSNP] Number of potential SNPs found: ",length(potentialSNPs))
		} else{
		    warning("[MethylToSNP] No potential SNPs found.")
		}
	}
	##########
	# STEP 2.
	#
	# Calculate SNP confidence
	# 
	
	snp.conf <- rep(NA, length(potentialSNPs))
	samples.low <- rep(NA, length(potentialSNPs))
	samples.mid <- rep(NA, length(potentialSNPs))
	samples.high <- rep(NA, length(potentialSNPs))
	
	for (i in 1:length(snp.conf)){
		count.low <- length(which(data[potentialSNPs[i],] <= 0.25))
		count.high <- length(which(data[potentialSNPs[i],] >= 0.75))
		count.mid <- length(data[potentialSNPs[i],]) - count.low - count.high

		L <- (count.high > 0) * (count.mid > 0) * (count.low > 0)
		snp.conf[i] <- round(L * (count.high + (count.mid / 2)) / dim(data)[2], 2)
		samples.low[i] <- count.low	
		samples.mid[i] <- count.mid
		samples.high[i] <- count.high
	}
	
	######
	# Concatenate with SNP information
	#
	# results <- NULL
	if (length(potentialSNPs) > 0) {
		results <- data.frame(
				row.names = potentialSNPs,
				confidence = snp.conf,
				samples_low = samples.low,
				samples_mid = samples.mid,
				samples_high = samples.high)		
		if (!missing(SNP)) {
			results <- cbind(results, SNP[potentialSNPs, ])
		}
		return(results)
	}
}


.checkMinfi <- function() {
	if (!exists('getBeta')) {
		if (!is_installed('minfi')) {
			stop("[MethylToSNP] Bioconductor minfi package is required to ")
		} else {
			library('minfi')
		}        	
	}
}

plotPotentialSNPs <- function(x, betas, horizontal=TRUE) {
	plotProbes(betas[rownames(x), ], horizontal)
}


plotProbes <- function(betas, horizontal=TRUE) {
	if (is(betas, "GenomicRatioSet") || is(betas, "GenomicMethylSet") || is(betas, "MethylSet") || is(betas, "RatioSet")) {
	    message("[MethylToSNP] Extracting beta values. Minfi required")
	    .checkMinfi()
	    betas <- getBeta(betas)
	}

	old_par <- par(no.readonly = TRUE)
	if (horizontal) {
		par(mar=c(6,5,2,2)) # bottom, left, top and right 
		stripchart(as.data.frame(t(betas)), ylab='Beta value', pch='-', vertical=TRUE, las=2, cex.axis=0.6) # las=horizontal)
	} else {
		par(mar=c(2,8,2,2))
		stripchart(as.data.frame(betas), xlab='Beta value', pch='|', las=2, cex=0.6)
	}
	par(old_par)
}

plotProbesMerged <- function(betas){
	old_par <- par(no.readonly = TRUE)
	par(mar=c(4,8,2,2))
	d <- dim(betas)
	plot(betas[1], ylim=c(0,1), xlab="Samples", ylab="Beta value")
	for (i in 2:d[1]) {
		plot(betas[i], add=TRUE)
	}
	par(old_par)
}

