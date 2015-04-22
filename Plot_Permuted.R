## Plot Heritability Estimates for GCTA ##
## January 27, 2015 ##
## Kristopher Standish ##

###################################################
## PARSE COMMAND LINE #############################
###################################################

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c( "/projects/janssen/Heritability/Manu_PhenoCovs_Derived.txt", "/projects/janssen/Heritability/20150316_Derived_MAF1_ALL_Manu_PhenoCovs_Derived", "TEMP", "10" )
Pheno_List <- LINE[1]
PathToData <- LINE[2]
Cohort_Name <- LINE[3]
Num_Perms <- as.numeric( LINE[4] )

print( "Running: Plot_Estimates.R" )
###################################################
## LOAD DATA ######################################
###################################################

## Load Phenotype List
print( "Loading Phenotype List")
PHENOS <- as.character( read.table( Pheno_List, header=F, fill=T )[,1] )

## Load Heritability Estimates for Phenotypes
print( "Loading/Compiling GCTA Results" )
EST <- VAR <- SE <- MOD <- list()
for ( p in 1:length(PHENOS) ) {
	pheno <- PHENOS[p]

	EST[[pheno]] <- list()
	VAR[[pheno]] <- array( , c(Num_Perms+1,4) ) ; colnames(VAR[[pheno]]) <- c("Vg","Ve","Vp","VgVp")
	SE[[pheno]] <- array( , c(Num_Perms+1,4) ) ; colnames(SE[[pheno]]) <- c("Vg","Ve","Vp","VgVp")
	MOD[[pheno]] <- array( , c(Num_Perms+1,6) ) ; colnames(MOD[[pheno]]) <- c("logL","logL0","LRT","df","Pval","n")
	rownames(VAR[[pheno]]) <- rownames(SE[[pheno]]) <- rownames(MOD[[pheno]]) <- c(paste("Perm",1:Num_Perms,sep="_"),"True")
	# Input Permuted Results
	for ( i in 1:Num_Perms ) {
		file_name <- paste(PathToData,"/4-PERM_",pheno,"_",i,".hsq",sep="")
		TEMP_TAB <- read.table( file_name, sep="\t",header=T,fill=T, colClasses=rep("character",3) )
		EST[[p]] <- TEMP_TAB
		## Pull out Variance Estimates for various parameters
		VAR[[pheno]][i,] <- as.numeric( TEMP_TAB[1:4,"Variance"] )
		SE[[pheno]][i,] <- as.numeric( TEMP_TAB[1:4,"SE"] )
		MOD[[pheno]][i,] <- as.numeric( TEMP_TAB[5:10,"Variance"] )
	}
	# Input True Results
	file_name <- paste(PathToData,"/3-REML_",pheno,".hsq",sep="")
	TEMP_TAB <- read.table( file_name, sep="\t",header=T,fill=T, colClasses=rep("character",3) )
	EST[[p]] <- TEMP_TAB
	## Pull out Variance Estimates for various parameters
	VAR[[pheno]][i+1,] <- as.numeric( TEMP_TAB[1:4,"Variance"] )
	SE[[pheno]][i+1,] <- as.numeric( TEMP_TAB[1:4,"SE"] )
	MOD[[pheno]][i+1,] <- as.numeric( TEMP_TAB[5:10,"Variance"] )

	
}
