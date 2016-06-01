# Permute Phenotype/Covarite Files #

# Permute_Phenos.R <Path/To/Pheno_Table> <Path/To/Cov_Table> <Perm_Number>

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c( paste("/projects/janssen/Heritability/20150316_Derived_MAF1_ALL_Manu_PhenoCovs_Derived/",c("PRE_MNw_DAS_FULL.txt","DEL_MNw_DAS_FULL.txt"),sep=""), 10 )
# LINE <- c( paste("/projects/janssen/Heritability/20160510_GCTA_Sing_MAF1_SNP_PC4_20160510_Single_Phenos/Phenos/",c("DEL_WAG20_DAS_FULL.txt","Cov_Ibl_DAS_PC4_FULL.txt"),sep=""), 10, "Ibl_DAS,PC1,PC2,PC3,PC4" )
Pheno_Table <- LINE[1]
Cov_Table <- LINE[2]
Perm_Number <- LINE[3]
Covs_List <- LINE[4]

## Print Inputs
print(paste( "Pheno_Table:", Pheno_Table ))
print(paste( "Cov_Table:", Cov_Table ))
print(paste( "Perm_Number:", Perm_Number ))
print(paste( "Covariates:", Covs_List ))

## Load Phenotype File
PHENO <- read.table(Pheno_Table,sep="\t",header=F)
print(paste( "dim(Pheno_Table):",dim(PHENO) ))

## Load Covariate File
if ( file.exists( Cov_Table ) ) { 
	COV <- read.table(Cov_Table,sep="\t",header=F)
	Covs_List.spl <- strsplit(Covs_List,",")[[1]]
	if ( any(grepl("^PC",Covs_List.spl)) ) {
		which_PCs <- 2 + grep("^PC",Covs_List.spl)
	}
	print(paste( "dim(Cov_Table):",dim(COV) ))
}

## Get Common Samples (just to be sure)
if ( file.exists( Cov_Table ) ) { 
	SAMPS <- intersect( as.character(PHENO[,1]), as.character(COV[,1]) )
	# PHENO <- PHENO[ which(PHENO[,1] %in% SAMPS), ]
	# COV <- COV[ which(COV[,1] %in% SAMPS), ]
	PHENO <- PHENO[ match(SAMPS,PHENO[,1]), ]
	COV <- COV[ match(SAMPS,COV[,1]), ]
}else{ SAMPS <- as.character(PHENO[,1]) }

## Permute ##
for ( p in 1:Perm_Number ) {
	## Re-Sample Patients
	reorder <- sample( 1:length(SAMPS) )
	REORDER <- SAMPS[ reorder ]

	## Compile/Write Permuted Phenotype Table
	 # Shuffle IDs, Leave Response Pheno
	PHENO.p <- PHENO
	PHENO.p[,1] <- PHENO.p[,2] <- REORDER
	 # Write Permuted Table
	Pheno_Out <- gsub( "txt", paste(p,"txt",sep="."), Pheno_Table )
	write.table( PHENO.p ,Pheno_Out,sep="\t",row.names=F,col.names=F,quote=F)
	
	## Compile/Write Permuted Covariate Table
	 # Shuffle IDs & PCs, Leave Other Covariates
	if ( file.exists( Cov_Table ) ) { 
		COV.p <- COV
		# Re-order IDs
		COV.p[,1] <- COV.p[,2] <- REORDER
		# Re-order PCs (if used)
		if ( exists("which_PCs") ) { COV.p[,which_PCs] <- COV.p[reorder,which_PCs] }
		# Write Permuted Table
		Cov_Out <- gsub( "txt", paste(p,"txt",sep="."), Cov_Table )
		write.table( COV.p ,Cov_Out,sep="\t",row.names=F,col.names=F,quote=F)
	}
}

######### END OF DOC ###########


