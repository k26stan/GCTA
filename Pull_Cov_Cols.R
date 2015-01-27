# Make Covariate File w/ EigenVectors #

# Make_Cov_Tab.R <Covar_1,Covar_2,Covar_3> <Path/To/EigenVectors> <Path/To/Covariate_Table> <Set>

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c(""/)
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
if ( grepl("PC",Covar_List) ) {
	USE_PCS <- 1
	PC_Count <- length( grep("PC",Covar_List) )
}else{ USE_PCS <- 0 }
print(paste( "USE_PCS:",USE_PCS ))

## Specify/Update Covariate List
if ( grepl(",",Covar_List) ) {
	Covar_List.2 <- sapply( strsplit( Covar_List, "," ), "[" )[,1]
}else{ Covar_List.2 <- Covar_List }
print( "Covar_List.2:" )
print( Covar_List.2 )

## Load PCs
if ( USE_PCS==1 ) {
	
	## Get rid of "PC" in Covar List if it's included
	Covar_List.2 <- Covar_List.2[-grep("PC",Covar_List.2)]
	print( "Covar_List.2:" )
	print( Covar_List.2 )

	## Load EigenVectors
	VEC <- read.table(PathToVec,header=F)
	print( "dim(EigenVectors):" )
	print( dim(VEC) )
}

## Merge Tables
COV_COLS <- c("IID",Covar_List.2)
COV.2 <- COV[,COV_COLS]
if ( USE_PCS==1 ) {
	VEC.2 <- VEC[,c(1,3:(PC_Count+3))]
	MRG <- merge(x=COV.2,y=VEC.2,by.x="IID",by.y="V1")
	print(paste( "dim(MRG):",dim(MRG) ))
}

## Write Table
if ( USE_PCS==1 ) { write.table(MRG[,c(1,1:ncol(MRG))],PathToNewCov,sep="\t",row.names=F,col.names=F,quote=F) }
if ( USE_PCS==0 ) { write.table(COV.2[,c(1,1:ncol(COV.2))],PathToNewCov,sep="\t",row.names=F,col.names=F,quote=F) }

######### END OF DOC ###########


