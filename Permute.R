# Permute Phenotype/Covarite Files #

# Permute_Phenos.R <Path/To/Pheno_Table> <Path/To/Cov_Table> <Perm_Number>

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c( paste("/projects/janssen/Heritability/20150316_Derived_MAF1_ALL_Manu_PhenoCovs_Derived/",c("PRE_MNw_DAS_FULL.txt","DEL_MNw_DAS_FULL.txt"),sep=""), 10 )
Pheno_Table <- LINE[1]
Cov_Table <- LINE[2]
Perm_Number <- LINE[3]

## Print Inputs
print(paste( "Pheno_Table:", Pheno_Table ))
print(paste( "Cov_Table:", Cov_Table ))
print(paste( "Perm_Number:", Perm_Number ))

## Load Phenotype File
PHENO <- read.table(Pheno_Table,sep="\t",header=F)
print(paste( "dim(Pheno_Table):",dim(PHENO) ))

## Load Covariate File
COV <- read.table(Cov_Table,sep="\t",header=F)
print(paste( "dim(Cov_Table):",dim(COV) ))

## Get Common Samples (just to be sure)
SAMPS <- intersect( as.character(PHENO[,1]), as.character(COV[,1]) )
PHENO <- PHENO[ which(PHENO[,1] %in% SAMPS), ]
COV <- COV[ which(COV[,1] %in% SAMPS), ]

## Permute ##
for ( p in 1:Perm_Number ) {
	## Sample Number of Patients
	REORDER <- sample( SAMPS )

	## Make New Tables
	PHENO.p <- PHENO
	PHENO.p[,1] <- PHENO.p[,2] <- REORDER
	COV.p <- COV
	COV.p[,1] <- COV.p[,2] <- REORDER

	## Path to Write
	Pheno_Out <- gsub( "txt", paste(p,"txt",sep="."), Pheno_Table )
	Cov_Out <- gsub( "txt", paste(p,"txt",sep="."), Cov_Table )

	## Write Table
	write.table( PHENO.p ,Pheno_Out,sep="\t",row.names=F,col.names=F,quote=F)
	write.table( COV.p ,Cov_Out,sep="\t",row.names=F,col.names=F,quote=F)

}

######### END OF DOC ###########


