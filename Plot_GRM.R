## Plot GRM Estimates for GCTA ##
## January 27, 2015 ##
## Kristopher Standish ##

###################################################
## PARSE COMMAND LINE #############################
###################################################

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- "/projects/janssen/Psych/GCTA/20150126b_PHENO_NAMES/1-GRM_FULL"
# LINE <- "/Users/kstandis/Downloads/1-GRM_FULL"
# LINE <- "/projects/janssen/Heritability/20151208_TestRun_Manu_PhenoCovs_Derived/1_GRM/1-GRM_FULL"
# LINE <- "/projects/janssen/Heritability/20151210_TestRun_Manu_PhenoCovs_Derived/1_GRM/1-GRM_FULL"
# LINE <- "/projects/janssen/Heritability/20151212_TestRun_LD.2_Manu_PhenoCovs_Derived/1_GRM/1-GRM_FULL"
PathToFile <- LINE[1]

###################################################
## LOAD DATA ######################################
###################################################

## Load Library to Plot
library(gplots) # heatmap.2

## Assign Function to Load Binary Genetic Relationship Files
ReadGRMBin <- function(prefix, AllN=F, size=4) {
	sum_i <- function(i) {
		return(sum(1:i))
	}
	BinFileName <- paste(prefix,".grm.bin",sep="")
	NFileName <- paste(prefix,".grm.N.bin",sep="")
	IDFileName <- paste(prefix,".grm.id",sep="")
	id <- read.table(IDFileName)
	n <- dim(id)[1]
	BinFile <- file(BinFileName, "rb")
	grm <- readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
	NFile <- file(NFileName, "rb")
	if(AllN==T){
		print( "AllN == T" )
		N <- readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
	}else{
		print( "AllN != T" )
		N <- readBin(NFile, n=1, what=numeric(0), size=size)
	}
	i <- sapply(1:n, sum_i)
	closeAllConnections()
	return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}

## Loop Through and Load Genetic Relationship Estimates
# PathToFile.CUT <- paste( PathToFile,"CUT",sep="." )
PathToFile.CUT <- gsub( "GRM_FULL","GRM_CUT",PathToFile )
AllN <- T
size <- 4
DAT.FULL <- ReadGRMBin( PathToFile, AllN, size )
DAT.CUT <- ReadGRMBin( PathToFile.CUT, AllN, size )

DAT.LIST <- list( DAT.FULL, DAT.CUT )
names(DAT.LIST) <- c("FULL","CUT")

#################################################
## REFORMAT & PLOT DATA #########################
#################################################

for ( d in 1:length(DAT.LIST) ) {
	DAT <- DAT.LIST[[d]]
	file_name <- names(DAT.LIST)[d]

	## Set Diagonal Matrix and then Fill
	DAT_ARR <- diag(DAT$diag) # * DAT$diag
	DAT_ARR[which(upper.tri(DAT_ARR))] <- DAT$off
	DAT_ARR[which(lower.tri(DAT_ARR))] <- t( DAT_ARR )[which(lower.tri(DAT_ARR))]
	rownames(DAT_ARR) <- colnames(DAT_ARR) <- DAT$id[,1]
	# X <- 1:10
	# DAT_ARR[X,X]

	## Specify plotting parameters
	COLORS <- c("firebrick2","sienna1","yellow1","black","chartreuse1","deepskyblue1","slateblue2")
	N_COLS <- 100
	COLS <- colorRampPalette(COLORS)(N_COLS)
	BRKS <- seq( -1,1,length.out=N_COLS+1)

	## Plot it
	PLOT_PRC <- 1
	X <- sample( 1:nrow(DAT_ARR), PLOT_PRC*nrow(DAT_ARR), replace=F )
	jpeg( paste(PathToFile,"_HEAT_",file_name,".jpeg",sep=""), height=2500, width=2500, pointsize=40 )
	heatmap.2(DAT_ARR[X,X], main=paste("Genetic Relationship Matrix -",100*PLOT_PRC,"% Cohort"), scale="none", symm=T, dendrogram="both", Rowv=T, Colv=T, col=COLS, breaks=BRKS, sepwidth=c(0,0), trace="none" )
	dev.off()

} # Close Data Loop

#################################################
## END OF DOC ###################################
#################################################


















