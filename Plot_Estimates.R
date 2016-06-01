## Plot Heritability Estimates for GCTA ##
## January 27, 2015 ##
## Kristopher Standish ##

###################################################
## PARSE COMMAND LINE #############################
###################################################

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c( "/projects/janssen/Pheno/PHENO_NAMES.txt", "/projects/janssen/Psych/GCTA/20150126b_PHENO_NAMES", "Cohort_Name" )
# LINE <- c( "/Users/kstandis/Downloads/20150126b_PHENO_NAMES/PHENO_NAMES.txt", "/Users/kstandis/Downloads/20150126b_PHENO_NAMES", "Cohort_Name" )
# LINE <- c( "/projects/janssen/Heritability/Manuscript_Phenos_Covs.txt", "/projects/janssen/Heritability/20150310_ALL2b_Manuscript_Phenos_Covs", "TEMP" )
# LINE <- c( "/projects/janssen/Heritability/Manuscript_Phenos_Covs.txt", "/projects/janssen/Heritability/20150310_ALL2b_Manuscript_Phenos_Covs", "TEMP" )
# LINE <- c( "/projects/janssen/Heritability/Manu_PhenoCovs_Single.SEL.txt", "/projects/janssen/Heritability/20150520_Single_SEL_MAF1_ALL_Manu_PhenoCovs_Single.SEL", "Full Cohort (MAF>1, SNP+Indel)" )
# LINE <- c( "/projects/janssen/Heritability/Manu_PhenoCovs_Derived.SEL.txt", "/projects/janssen/Heritability/20150520_Derived_SEL_MAF1_ALL_Manu_PhenoCovs_Derived.SEL", "Full Cohort (MAF>1, SNP+Indel)" )
# LINE <- c( "/projects/janssen/Heritability/Manu_PhenoCovs_Derived.txt", "/projects/janssen/Heritability/20160125_GCTA_Derived_MAF1_SNP_Manu_PhenoCovs_Derived", "Full Cohort (MAF>1, SNP)" )
Pheno_List <- LINE[1]
PathToData <- LINE[2]
Cohort_Name <- LINE[3]

## Print Inputs
print( "!!Running: Plot_Estimates.R" )
print(paste( "Pheno List:", Pheno_List ))
print(paste( "Path:", PathToData ))
print(paste( "Cohort Tag:", Cohort_Name ))

###################################################
## LOAD DATA ######################################
###################################################

## Load Phenotype List
print( "Loading Phenotype List")
PHENOS <- as.character( read.table( Pheno_List, header=F, fill=T )[,1] )
PHENOS.n <- PHENOS

## ***************** OPTIONAL: SINGLE ************************* ##

# ## Change Names/Order for SINGLE Timepoint Estimates
# PHENOS.n[1:28] <- paste( "POST",PHENOS.n[1:28],sep="_" )
# ReOrder <- seq(1,28,7) + rep( 1:7-1,rep(4,7) )
# PHENOS.n <- PHENOS.n[ c(ReOrder,29:length(PHENOS.n)) ]
# PHENOS <- PHENOS[ c(ReOrder,29:length(PHENOS)) ]

# ## To Plot only Subset of Phenotypes
#  # "POST" treatment Measurements
# PHENOS <- PHENOS[ grep("POST",PHENOS.n) ]
# PHENOS.n <- PHENOS.n[ grep("POST",PHENOS.n) ]
#  # "DEL" Measurements
# PHENOS <- PHENOS[ grep("DEL",PHENOS.n) ]
# PHENOS.n <- PHENOS.n[ grep("DEL",PHENOS.n) ]

## ***************** OPTIONAL: DERIVED ************************* ##

# ## To Plot only Subset of Phenotypes
 # Remove "28"JC Phenos
PHENOS <- PHENOS[ grep("JC28",PHENOS, invert=T) ]
PHENOS.n <- PHENOS.n[ grep("JC28",PHENOS.n, invert=T) ]
 # Extract only "DEL" Phenos
PHENOS <- PHENOS[ grep("DEL",PHENOS) ]
PHENOS.n <- PHENOS.n[ grep("DEL",PHENOS.n) ]
#  # Extract only "MNa" Phenos
# PHENOS <- PHENOS[ grep("MNa",PHENOS) ]
# PHENOS.n <- PHENOS.n[ grep("MNa",PHENOS.n) ]
#  # Extract non-"MNa" Phenos
# PHENOS <- PHENOS[ grep("MNa",PHENOS, invert=T) ]
# PHENOS.n <- PHENOS.n[ grep("MNa",PHENOS.n, invert=T) ]



## Load Heritability Estimates for Phenotypes
print( "Loading/Compiling GCTA Results" )
EST <- list()
VAR <- array( , c(length(PHENOS),4) ) ; colnames(VAR) <- c("Vg","Ve","Vp","VgVp")
SE <- array( , c(length(PHENOS),4) ) ; colnames(SE) <- c("Vg","Ve","Vp","VgVp")
MOD <- array( , c(length(PHENOS),6) ) ; colnames(MOD) <- c("logL","logL0","LRT","df","Pval","n")
rownames(VAR) <- rownames(SE) <- rownames(MOD) <- gsub(".txt","",PHENOS)
for ( p in 1:length(PHENOS) ) {
	pheno <- PHENOS[p]
	file_name <- paste(PathToData,"/3_REML/3-REML_",pheno,".hsq",sep="")
	TEMP_TAB <- read.table( file_name, sep="\t",header=T,fill=T, colClasses=rep("character",3) )
	EST[[p]] <- TEMP_TAB
	## Pull out Variance Estimates for various parameters
	VAR[p,] <- as.numeric( TEMP_TAB[1:4,"Variance"] )
	SE[p,] <- as.numeric( TEMP_TAB[1:4,"SE"] )
	MOD[p,] <- as.numeric( TEMP_TAB[5:10,"Variance"] )
}

###################################################
## PLOT VARIANCE ESTIMATES ########################
###################################################
print( "Plotting Results" )

## Basic Plot for (all) Phenotypes
COLS.list <- c("firebrick1","gold1","chartreuse1","dodgerblue1") # c("firebrick2","gold2","chartreuse2","deepskyblue2","slateblue2")
COLS <- rep( COLS.list[1], nrow(VAR) )
COLS[ grep("SJC",rownames(VAR)) ] <- "gold4"
COLS[ grep("TJC",rownames(VAR)) ] <- "dodgerblue4"
COLS[ grep("CRP",rownames(VAR)) ] <- "chartreuse4"
COLS[ grep("rSJC",rownames(VAR)) ] <- "gold1"
COLS[ grep("rTJC",rownames(VAR)) ] <- "dodgerblue1"
COLS[ grep("lCRP",rownames(VAR)) ] <- "chartreuse1"
LTYS <- rep( 2, nrow(VAR) )
LTYS[ grep("DEL",rownames(VAR)) ] <- 1
PCHS <- rep( 1, nrow(VAR) )
PCHS[ grep("DEL",rownames(VAR)) ] <- 20

## Set Plot Parameters
YLIM <- c(min( 0,VAR[,"VgVp"]-SE[,"VgVp"], na.rm=T), max(1,max(VAR[,"VgVp"]+SE[,"VgVp"],na.rm=T)) )
XLIM <- c( 0,nrow(VAR)+1 )
WHICH_SIG <- which( MOD[,"Pval"] < .05 )
## Open Plot
jpeg( paste(PathToData,"/GCTA_Estimates.ALL.jpeg",sep=""), height=1500,width=1000+50*length(PHENOS), pointsize=30 )
plot( 0,0,type="n", xlim=XLIM, ylim=YLIM+c(0,.4), main=paste("Heritability Estimate -",Cohort_Name), ylab="% Phenotypic Variance", xlab="Phenotype", yaxt="n", xaxt="n" )
 # Vertical Grid Lines
abline( h=seq(-2,XLIM[2]+.4,.1), lty=2, col="grey50", lwd=1 )
abline( h=c(0,1), lty=1, col="black", lwd=1 )
 # Horizontal Lines/Data
arrows( 1:nrow(VAR), 0, 1:nrow(VAR), 1, col="black", lwd=1, length=0 )
arrows( 1:nrow(VAR), VAR[,"VgVp"]-SE[,"VgVp"], 1:nrow(VAR), VAR[,"VgVp"]+SE[,"VgVp"], col=COLS, lty=LTYS, code=3, angle=90, length=.15, lwd=6 )
points( 1:nrow(VAR), VAR[,"VgVp"], col=COLS, pch=PCHS, cex=1.4, lwd=3 )
 # Axis/Labels/Significance
axis(2, at=seq(-2,2,.2) )
text( 1:nrow(VAR)-.25, .05+sapply( VAR[,"VgVp"]+SE[,"VgVp"], function(x) max(x,1) ), labels=PHENOS.n, pos=4, cex=.8, col=COLS, srt=90 )
if ( length(WHICH_SIG)>0 ) { text( WHICH_SIG, rep(-.05,length(WHICH_SIG)), labels="*", col=COLS[WHICH_SIG], cex=1.5 ) }
dev.off()

###################################################
## SAVE COMPILED TABLES ###########################
###################################################

## Compile Data
COMPILE <- list(EST=EST,VAR=VAR,SE=SE,MOD=MOD)
PathToSave <- paste(PathToData,"/GCTA_Estimates.ALL.Rdata",sep="")
save( COMPILE, file=PathToSave )

###################################################
## END OF DOC #####################################
###################################################






# ## Set Plot Parameters
# XLIM <- c(min( 0,VAR[,"VgVp"]-SE[,"VgVp"], na.rm=T), max(1,max(VAR[,"VgVp"]+SE[,"VgVp"],na.rm=T)) )
# YLIM <- c( 0,nrow(VAR)+1 )
# WHICH_SIG <- which( MOD[,"Pval"] < .05 )
# ## Open Plot
# jpeg( paste(PathToData,"/GCTA_Estimates.jpeg",sep=""), width=1500,height=1000+50*length(PHENOS), pointsize=30 )
# plot( 0,0,type="n", xlim=XLIM+c(0,.4), ylim=YLIM, main=paste("Heritability Estimate -",Cohort_Name), xlab="% Phenotypic Variance", ylab="Phenotype", yaxt="n", xaxt="n" )
#  # Vertical Grid Lines
# abline( v=seq(-2,XLIM[2]+.4,.1), lty=2, col="grey50", lwd=1 )
# abline( v=c(0,1), lty=1, col="black", lwd=1 )
#  # Horizontal Lines/Data
# arrows( 0, 1:nrow(VAR), 1, 1:nrow(VAR), col="black", lwd=1, length=0 )
# arrows( VAR[,"VgVp"]-SE[,"VgVp"], 1:nrow(VAR), VAR[,"VgVp"]+SE[,"VgVp"], 1:nrow(VAR), col=COLS, lty=LTYS, code=3, angle=90, length=.15, lwd=6 )
# points( VAR[,"VgVp"], 1:nrow(VAR), col=COLS, pch=PCHS, cex=1.4, lwd=3 )
#  # Axis/Labels/Significance
# axis(1, at=seq(-2,2,.2) )
# text( sapply( VAR[,"VgVp"]+SE[,"VgVp"], function(x) max(x,1) ) , 1:nrow(VAR), labels=gsub("LT8_DEL_","",rownames(VAR)), pos=4, cex=.8, col=COLS )
# if ( length(WHICH_SIG)>0 ) { text( rep(-.05,length(WHICH_SIG)), WHICH_SIG, labels="*", col=COLS[WHICH_SIG], cex=1.5 ) }
# dev.off()
