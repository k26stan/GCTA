## Plot Heritability Estimates for GCTA ##
## January 27, 2015 ##
## Kristopher Standish ##

###################################################
## PARSE COMMAND LINE #############################
###################################################

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c( "/projects/janssen/Pheno/PHENO_NAMES.txt", "/projects/janssen/Psych/GCTA/20150126b_PHENO_NAMES", "Cohort_Name" )
# LINE <- c( "/Users/kstandis/Downloads/20150126b_PHENO_NAMES/PHENO_NAMES.txt", "/Users/kstandis/Downloads/20150126b_PHENO_NAMES", "Cohort_Name" )
Pheno_List <- LINE[1]
PathToData <- LINE[2]
Cohort_Name <- LINE[3]

print( "Running: Plot_Estimates.R" )
###################################################
## LOAD DATA ######################################
###################################################

## Load Phenotype List
print( "Loading Phenotype List")
PHENOS <- as.character( read.table( Pheno_List, header=F )[,1] )

## Load Heritability Estimates for Phenotypes
print( "Loading/Compiling GCTA Results" )
EST <- list()
VAR <- array( , c(length(PHENOS),4) ) ; colnames(VAR) <- c("Vg","Ve","Vp","VgVp")
SE <- array( , c(length(PHENOS),4) ) ; colnames(SE) <- c("Vg","Ve","Vp","VgVp")
MOD <- array( , c(length(PHENOS),6) ) ; colnames(MOD) <- c("logL","logL0","LRT","df","Pval","n")
rownames(VAR) <- rownames(SE) <- rownames(MOD) <- gsub(".txt","",PHENOS)
for ( p in 1:length(PHENOS) ) {
	pheno <- PHENOS[p]
	file_name <- paste(PathToData,"3-REML_",pheno,".hsq",sep="")
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
COLS <- "steelblue2" # c("firebrick2","gold2","chartreuse2","deepskyblue2","slateblue2")

## Set Plot Parameters
XLIM <- c(min( VAR[,"VgVp"]-SE[,"VgVp"], na.rm=T), max(1,max(VAR[,"VgVp"]+SE[,"VgVp"],na.rm=T)) )
YLIM <- c( 0,nrow(VAR)+1 )
WHICH_SIG <- which( MOD[,"Pval"] < .05 )
## Open Plot
jpeg( paste(PathToData,"/GCTA_Estimates.jpeg",sep=""), width=1500,height=1000+50*length(PHENOS), pointsize=30 )
plot( 0,0,type="n", xlim=XLIM+c(0,.4), ylim=YLIM, main=paste("Heritability Estimate -",Cohort_Name), xlab="% Phenotypic Variance", ylab="Phenotype", yaxt="n", xaxt="n" )
 # Vertical Grid Lines
abline( v=seq(-2,XLIM[2]+.4,.1), lty=2, col="grey50", lwd=1 )
abline( v=c(0,1), lty=1, col="black", lwd=1 )
 # Horizontal Lines/Data
arrows( 0, 1:nrow(VAR), 1, 1:nrow(VAR), col="black", lwd=1, length=0 )
arrows( VAR[,"VgVp"]-SE[,"VgVp"], 1:nrow(VAR), VAR[,"VgVp"]+SE[,"VgVp"], 1:nrow(VAR), col=COLS, code=3, angle=90, length=.15, lwd=6 )
points( VAR[,"VgVp"], 1:nrow(VAR), col=COLS, pch=20, cex=1.4 )
axis(1, at=seq(-2,2,.2) )
 # 
text( sapply( VAR[,"VgVp"]+SE[,"VgVp"], function(x) max(x,1) ) , 1:nrow(VAR), labels=gsub("LT8_DEL_","",rownames(VAR)), pos=4, cex=.8, col=COLS )
if ( length(WHICH_SIG)>0 ) { text( rep(-.05,length(WHICH_SIG)), WHICH_SIG, labels="*", col=COLS, cex=1.5 ) }
dev.off()

###################################################
## END OF DOC #####################################
###################################################
