# Make Covariate File w/ EigenVectors #

# Make_Cov_Tab.R <Covar_1,Covar_2,Covar_3> <Path/To/EigenVectors> <Path/To/Covariate_Table> <Set>

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("DEL_MNa_DAS,PRE_MNa_DAS,PC1,PC2,PC3","/projects/janssen/Heritability/20151208_TestRun_Manu_PhenoCovs_Derived/2-PCA.eigenvec","/projects/janssen/ASSOCIATION/PH-PHENOTYPES/20150506_Derived_Pheno_Table.txt","/projects/janssen/Heritability/20151208_TestRun_Manu_PhenoCovs_Derived/Phenos/DEL_MNa_DAS PRE_MNa_DAS_PC1_PC2_PC3_FULL.txt")
# LINE <- c(",PC1,PC2,PC3","/projects/janssen/Heritability/20151208_TestRun_Manu_PhenoCovs_Derived/2-PCA.eigenvec","/projects/janssen/ASSOCIATION/PH-PHENOTYPES/20150506_Derived_Pheno_Table.txt","/projects/janssen/Heritability/20151208_TestRun_Manu_PhenoCovs_Derived/Phenos/_PC1_PC2_PC3_FULL.txt")
Covar_List <- LINE[1]
PathToVec <- LINE[2]
PathToCov <- LINE[3]
PathToNewCov <- LINE[4]

## Print Inputs
print(paste( "Covar_List:", Covar_List ))
print(paste( "PathToVec:", PathToVec ))
print(paste( "PathToCov:", PathToCov ))
print(paste( "PathToNewCov:", PathToNewCov ))

## Load Covariate File
COV <- read.table(PathToCov,sep="\t",header=T)
print(paste( "dim(CovarTable):",dim(COV) ))

## Use PCs?
if ( any(grepl("PC",Covar_List)) ) {
	USE_PCS <- 1
}else{ USE_PCS <- 0 }
print(paste( "USE_PCS:",USE_PCS ))

## Specify/Update Covariate List
if ( grepl(",",Covar_List) ) {
	Covar_List.2 <- strsplit( Covar_List, "," )[[1]]
}else{ Covar_List.2 <- Covar_List }
print( "Covar_List.2:" )
print( Covar_List.2 )

## Load PCs
if ( USE_PCS==1 ) {
	# Get PC Count
	PC_Count <- length( grep("PC",Covar_List.2) )
	print(paste( "PC_Count:",PC_Count ))
	# Get rid of "PC" in Covar List if it's included
	Covar_List.2 <- Covar_List.2[-grep("PC",Covar_List.2)]
	print( "Covar_List.2:" )
	print( Covar_List.2 )
	# Load EigenVectors
	VEC <- read.table(PathToVec,header=F)
	print( "dim(EigenVectors):" )
	print( dim(VEC) )
}

## Check for Remaining Covariates
if ( Covar_List.2[1]=="" ) {
	OUT <- VEC[,c(1,3:(PC_Count+2))]
}else{
	# Pull out Covariate Columns
	COV_COLS <- c("IID",Covar_List.2)
	COV.2 <- COV[,COV_COLS]
	# Merge w/ PCs
	if ( USE_PCS==1 ) {
		VEC.2 <- VEC[,c(1,3:(PC_Count+2))]
		OUT <- merge(x=COV.2,y=VEC.2,by.x="IID",by.y="V1")
	}else{
		OUT <- COVS.2
	}
}
print(paste( "dim(OUT):",dim(OUT) ))

## Write Table
write.table(OUT[,c(1,1:ncol(OUT))],PathToNewCov,sep="\t",row.names=F,col.names=F,quote=F)
# if ( USE_PCS==1 ) { write.table(OUT[,c(1,1:ncol(OUT))],PathToNewCov,sep="\t",row.names=F,col.names=F,quote=F) }
# if ( USE_PCS==0 ) { write.table(COV.2[,c(1,1:ncol(COV.2))],PathToNewCov,sep="\t",row.names=F,col.names=F,quote=F) }

######### END OF DOC ###########


