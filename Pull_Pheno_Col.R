# Make Phenotype File #

# Make_Pheno_Tab.R <Path/To/Full_Table> <Pheno_Name> <Path/To/New_Pheno_Table>

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c("/projects/janssen/ASSOCIATION/PH-PHENOTYPES/20150506_Derived_Pheno_Table.txt","DEL_MNa_DAS","/projects/janssen/Heritability/20151208_TestRun_Manu_PhenoCovs_Derived/Phenos/DEL_MNa_DAS_FULL.txt")
Pheno_Table <- LINE[1]
Pheno_Name <- LINE[2]
New_Pheno_Table <- LINE[3]

## Print Inputs
print( "!!Running: Pull_Pheno_Cols.R" )
print(paste( "Pheno_Table:", Pheno_Table ))
print(paste( "Pheno_Name:", Pheno_Name ))
print(paste( "New_Pheno_Table:", New_Pheno_Table ))

## Load Phenotype File
PHENO <- read.table(Pheno_Table,sep="\t",header=T)
# print(paste( "dim(Pheno_Table):",dim(PHENO) ))

## Pull Phenotype Columne
PHENO_COL <- c("IID",Pheno_Name)
PHENO.2 <- PHENO[,PHENO_COL]

## Remove Missing Values
PHENO.3 <- PHENO.2[ which(!is.na(PHENO.2[,Pheno_Name])) ,]

## Write Table
write.table(PHENO.3[,c(1,1:ncol(PHENO.3))],New_Pheno_Table,sep="\t",row.names=F,col.names=F,quote=F)

######### END OF DOC ###########


