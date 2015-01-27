# Make Phenotype File #

# Make_Pheno_Tab.R <Path/To/Full_Table> <Pheno_Name> <Path/To/New_Pheno_Table>

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c(""/)
Pheno_Table <- LINE[1]
Pheno_Name <- LINE[2]
New_Pheno_Table <- LINE[3]

## Print Inputs
print(paste( "Pheno_Table:", Pheno_Table ))
print(paste( "Pheno_Name:", Pheno_Name ))
print(paste( "New_Pheno_Table:", New_Pheno_Table ))

## Load Covariate File
PHENO <- read.table(Pheno_Table,sep="\t",header=T)
print(paste( "dim(Pheno_Table):",dim(PHENO) ))

## Merge Tables
PHENO_COL <- c("IID",Pheno_Name)
PHENO.2 <- PHENO[,PHENO_COL]

## Write Table
write.table(PHENO.2[,c(1,1:ncol(PHENO.2))],New_Pheno_Table,sep="\t",row.names=F,col.names=F,quote=F)

######### END OF DOC ###########


