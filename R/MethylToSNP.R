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
#' @param SNP SNP annotation
#' @param verbose show additional information. Useful for debugging.
#' @return Detected probes with 3-tier SNP-like methylation pattern along with their reliability scores and SNP annotation
#' @examples
#' MethylToSNP(data)
#' @export
MethylToSNP <- function(data, gap.ratio = 0.75, gap.sum.ratio = 0.5, verbose=FALSE, outlier.sd = 3.0, SNP
	method="clustering")
{
    if ( (!length(dim(data)) == 2) || (dim(data)[2] < 2)){
        stop("[MethylToSNP] There have to be at least 2 samples")
    }

    if (dim(data)[2] < 50){
        message("[MethylToSNP] Warning, SNP detection may be unreliable in datasets with less than 50 samples")
    }
    
    if (method != "feature" && method != "clustering" ){
		stop("[MethylToSNP] You have entered an unacceptable method name. Please use only 'feature' or 'clustering'.")
    }
    
	if (gap.ratio >= 1 || gap.ratio <= 0){
		stop("[MethylToSNP] You have entered an unacceptable gap.ratio value. This must be a number between zero and one.")
	}
	
    if (gap.sum.ratio <= 0 || gap.sum.ratio >= 1){
		stop("[MethylToSNP] You have entered an unacceptable gap.sum.ratio value. This must be a number between zero and one.")
	}
    
    if (is(data, "GenomicRatioSet") || is(data, "GenomicMethylSet") || is(data, "MethylSet") || is(data, "RatioSet")) {
        # .isMatrixBackedOrStop(data, "MethylToSNP")
        if (verbose) 
            message("[MethylToSNP] Calculating beta matrix.")
        data <- getBeta(data)
    }
    
    # if (sum((beta >= 1) || (beta <= 1)) > 0) {
    #    warning('Detected incorrect beta values. Assuming input was M values and transforming to Beta values')
    #    beta <- (2^beta)/((2^beta)+1)
    # }
    
    ##########
    # STEP 1.
    # Identify probes that could be potential SNPs in a given set of samples
    
	potentialSNPs <- c()
	
	# for each probe f
	for (i in 1: dim(data)[1]){
		data.range <- diff(range(data[i,], na.rm = TRUE))
		if (abs(data.range) >= 0.5){
			
			if (method == "feature"){	
				pofi <- as.numeric(sort(data[i,]))
				p.len <- length(pofi)
				sort.diff <- rep(NA, (p.len-1))
			
				for(j in 1:(p.len - 1)){
					sort.diff[j] <- pofi[j+1] - pofi[j]
				}
				
				sort.sort.diff <- sort(sort.diff, decreasing = 	TRUE)[1:2]
				
			} else if(method == "clustering"){
				# if(sum(is.na(data[i,])) != 0){
				# 	pofi <- data[i,-which(is.na(data[i,]))]
				# } else{
				# 	pofi <- data[i,]
				# }
				data_for_clustering <- t(data[i,])

				kofp <- Ckmeans.1d.dp(x=data_for_clustering, y=weights, k=3)


				# kofp <- Ckmedian.1d.dp(x=data_for_clustering, k=3)
				# kofp <- Ckmedian.1d.dp(x=t(pofi), k=3)
				# print(kofp)
				# plot(kofp)

				top <- which(kofp$centers == max(kofp$centers))
				bottom <- which(kofp$centers == min(kofp$centers))
				middle <- c(1:3)[-c(top, bottom)]

				if (outlier.sd != FALSE) {
					top_outliers <- which(abs(scale(data[i, kofp$cluster == top])) > outlier.sd)
					kofp$cluster[kofp$cluster == top][] <- 4  # assign outliers to non-existing cluster 4

					bottom_outliers <- which(abs(scale(data[i, kofp$cluster == bottom])) > outlier.sd)
					kofp$cluster[kofp$cluster == bottom][] <- 4  # assign outliers to non-existing cluster 4

					middle_outliers <- which(abs(scale(data[i, kofp$cluster == middle])) > outlier.sd)
					kofp$cluster[kofp$cluster == middle][] <- 4  # assign outliers to non-existing cluster 4
				}
				
				sort.diff <- min(data[i, which( kofp$cluster == top)], na.rm = TRUE) - max(data[i, which(kofp$cluster == middle)], na.rm = TRUE)
				sort.diff <- c(sort.diff, min(data[i, which( kofp$cluster == middle)], na.rm = TRUE) - max(data[i, which(kofp$cluster == bottom)], na.rm = TRUE))
				
				sort.sort.diff <- sort(sort.diff, decreasing = TRUE)
			}
					
			###
		    # Apply gap thresholds to decide whether to add a cg probe to the list of potential SNPs
		    #
			if ((sort.sort.diff[2] >= (gap.ratio*(sort.sort.diff[1]))) && (sum(sort.sort.diff) >= gap.sum.ratio*data.range)){
			    if (verbose) {
			        message(rownames(data)[i])
			        message(sort.sort.diff)
			    }
			    potentialSNPs <- append(potentialSNPs, rownames(data)[i])
			}
		}		
		
		###
		# Progress indicator
		#
		if (verbose){
			if((i %% 10000) == 0){
				message('[MethylToSNP] Processed', " ", i, " ", length(potentialSNPs))
			}
		}
	}
	###
	# Summary
	#
	if(verbose){
		message("[MethylToSNP] Finished SNP search: ", Sys.time())
		if (length(potentialSNPs) > 0 ){
		    message("[MethylToSNP] Number of potential SNPs found: ",length(potentialSNPs))
		} else{
		    message("[MethylToSNP] No potential SNPs found.")
		}
	}
	##########
	# STEP 2.
	# Now that we have a list of potential SNPs, which can check which ones of those are already known, and calculate a minimum true positive rate
	#
	
	# FIXME: THERE ARE NO KNOWN RS IN minfi !!!!
	# if there are we should remove them!

	# some of the probes on the chip are actually named as the SNPs they hit
	#knownrs <- grep("rs", rownames(data))
	
	# plenty of probes on the chip overlap other previously identified SNPs, directions for getting the most up-to-date version of this file are in the header information
	# known.list <- snp.list
	
	knownSNPslist <- 0
	potentials.noknown <- potentialSNPs
	
	# if (any(potentialSNPs %in% rownames(data)[knownrs])){
	# 	knownSNPslist <- potentialSNPs[potentialSNPs %in% rownames(data)[knownrs]]
	# }
	
	# if (known.list && any(potentialSNPs %in% known.list[,4])){
	# 	if (knownSNPslist[1] == 0){
	# 		knownSNPslist <- potentialSNPs[potentialSNPs %in% known.list[,4]]
	# 	} else{
	# 		knownSNPslist <- c(knownSNPslist, potentialSNPs[(potentialSNPs %in% known.list[,4])])
	# 	}
	# }
	
	# potentials.noknown <- potentials.noknown[-which(potentialSNPs %in% knownSNPslist)]
	
	# if (is.character(knownSNPslist[1])){
	# 	true.positive <- length(knownSNPslist) / length(potentialSNPs)
	# }else{
	# 	true.positive <- 0
	# }
	
	######
	# Calculate SNP confidence
	#
	snp.conf <- rep(NA, length(potentialSNPs))
	
	for (i in 1:length(snp.conf)){
		count.low <- length(which(data[potentialSNPs[i],] <= 0.25))
		count.high <- length(which(data[potentialSNPs[i],] >= 0.75))
		count.mid <- length(data[potentialSNPs[i],]) - count.low - count.high

		L <- (count.high == 0) * (count.mid == 0) * (count.low == 0)
		snp.conf[i] <- L * (count.high + (count.mid / 2)) / dim(data)[2]
		
		# if (count.low == 0 && count.mid == 0){
		# 	levels <- 0
		# } else if (count.low == 0 && count.high == 0){
		# 	levels <- 0
		# } else if (count.mid == 0 && count.high == 0){
		# 	levels <- 0
		# } else{
		#     levels <- 1
		# }
		# snp.conf[i] <- ((count.high + (count.mid / 2)) / dim(data)[2]) * levels
	}
	
	######
	# Check if SNPs in probe sequence
	#
	# TODO: That's an inefficient search, need to use genomic ranges
	# TODO: initialize with empty list and then append
	probes.w.SNPs <- 0
	for (i in 1:length(potentials.noknown)){
		# WHAT TO DO WITH ANNOTATIONS
		# if (annotations[potentials.noknown[i], which.anno.isSNPs] != ""){
		# 	if (probes.w.SNPs == 0){
		# 		probes.w.SNPs <- potentials.noknown[i]
		# 	} else {
		# 	    probes.w.SNPs <- c(probes.w.SNPs, potentials.noknown[i])
		#     }
		# }
	}

	all.res.out <- list(
		potentialSNPs = potentialSNPs,
		potentials.noknown = potentials.noknown,
		knownSNPs = knownSNPslist,
		probeswithSNPs = probes.w.SNPs,
		true.positive = true.positive,
		snp.conf = snp.conf)
	
	if (verbose){
		if (length(potentials.noknown) > 1){
		    message("[MethylToSNP] Finished checking reference SNP list. Number of potential SNPs not previously reported: ", length(potentials.noknown))
		} else if(length(potentials.noknown) == 1 && potentials.noknown[1] != 0){
		    message("[MethylToSNP] Finished checking reference SNP list. Number of potential SNPs not previously reported: ", length(potentials.noknown))
		} else{
		    message("[MethylToSNP] No potential new SNPs found.")
		}
	}
	
	return(invisible(all.res.out))
}




plotProbes <- function(betas, horizontal=TRUE){
	old_par <- par(no.readonly = TRUE)
	if (horizontal) {
		par(mar=c(6,5,2,2)) # bottom, left, top and right 
		stripchart(as.data.frame(t(betas)), ylab='Beta value', pch='-', vertical=TRUE, las=2, cex.axis=0.6, cex.names=0.8) # las=horizontal)
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



