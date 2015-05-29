## Plot Heritability Estimates for GCTA ##
## January 27, 2015 ##
## Kristopher Standish ##

###################################################
## PARSE COMMAND LINE #############################
###################################################

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- c( "/projects/janssen/Heritability/Manu_PhenoCovs_Derived.txt", "/projects/janssen/Heritability/20150316_Derived_MAF1_ALL_Manu_PhenoCovs_Derived", "TEMP", "1000" )
# LINE <- c( "/projects/janssen/Heritability/Manu_PhenoCovs_Derived.txt", "/projects/janssen/Heritability/20150316_Derived_MAF1_ALL_Manu_PhenoCovs_Derived", "1000" )
# LINE <- c( "/projects/janssen/Heritability/Manu_PhenoCovs_Derived.txt", "/projects/janssen/Heritability/20150506_Derived_MAF1_ALL_Manu_PhenoCovs_Derived", "100" )
# LINE <- c( "/projects/janssen/Heritability/Manu_PhenoCovs_Derived.SEL.txt", "/projects/janssen/Heritability/20150506_Derived_MAF1_ALL_Manu_PhenoCovs_Derived.SEL", "1000" )
# LINE <- c( "/projects/janssen/Heritability/Manu_PhenoCovs_Derived.SEL.txt", "/projects/janssen/Heritability/20150520_Derived_SEL_MAF1_ALL_Manu_PhenoCovs_Derived.SEL", "1000" )
# LINE <- c( "/projects/janssen/Heritability/Manu_PhenoCovs_Single.SEL.txt", "/projects/janssen/Heritability/20150520_Single_SEL_MAF1_ALL_Manu_PhenoCovs_Single.SEL", "1000" )
Pheno_List <- LINE[1]
PathToData <- LINE[2]
Num_Perms <- as.numeric( LINE[3] )

print( "Running: Plot_Estimates.R" )
###################################################
## LOAD DATA ######################################
###################################################

## Load Phenotype List
print( "Loading Phenotype List")
PHENOS.l <- read.table( Pheno_List, header=F, fill=T )
PHENOS <- as.character( PHENOS.l[,1] )
# PHENOS <- as.character( PHENOS.l[ which(PHENOS.l[,2]!=""),1] )
# TEMP.PHENO.list <- as.character( read.table( "TEMP.PHENO.list", header=F, fill=T )[,1] )
# TEMP.PHENO.list <- gsub("_1.hsq","",TEMP.PHENO.list)
# PHENOS <- PHENOS[ -grep("MNw",PHENOS) ]
# PHENOS <- PHENOS[ -grep("28",PHENOS) ]
# PHENOS.post <- PHENOS[ grep("POST",PHENOS) ]
# PHENOS.del <- PHENOS[ grep("DEL",PHENOS) ]

## ***************** OPTIONAL: SINGLE ************************* ##

## Change Names/Order for SINGLE Timepoint Estimates
ReOrder <- seq(1,28,7) + rep( 1:7-1,rep(4,7) )
PHENOS <- PHENOS[ c(ReOrder,29:length(PHENOS)) ]

## To Plot only Subset of Phenotypes
#  # "POST" treatment Measurements
# PHENOS <- PHENOS[ grep("POST",PHENOS) ]
 # "DEL" Measurements
PHENOS <- PHENOS[ grep("DEL",PHENOS) ]

## ***************** OPTIONAL: DERIVED ************************* ##

# ## To Plot only Subset of Phenotypes
#  # Remove "28"JC Phenos
# PHENOS <- PHENOS[ grep("28",PHENOS, invert=T) ]
#  # Extract only "DEL" Phenos
# PHENOS <- PHENOS[ grep("DEL",PHENOS) ]
#  # Extract only "MNa" Phenos
# PHENOS <- PHENOS[ grep("MNa",PHENOS) ]
#  # Extract non-"MNa" Phenos
# PHENOS <- PHENOS[ grep("MNa",PHENOS, invert=T) ]



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
		file_name <- paste(PathToData,"/4-PERM/",pheno,"_",i,".hsq",sep="")
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

	## Update Output
	print(paste( "Done with",pheno,":",p,"of",length(PHENOS) ))
}

## Compile and Save Data
COMPILE.full <- list( EST, VAR, SE, MOD )
names(COMPILE.full) <- c("EST","VAR","SE","MOD")
save( COMPILE.full, file=paste(PathToData,"/4-PERM_Compile.Rdata",sep="") )
# load( paste(PathToData,"/4-PERM_Compile.Rdata",sep="") )
COMPILE <- COMPILE.full
 # Optionally Remove 
for ( i in 1:length(COMPILE) ) {
	COMPILE[[i]] <- COMPILE[[i]][PHENOS]
}

###################################################
## MAKE PLOTS #####################################
###################################################
# PHENOS <- PHENOS.post
# COMPILE.post <- lapply( COMPILE.full, function(x) which(x %in% PHENOS.post) )

## Basic Plot for (all) Phenotypes
COLS.list <- c("firebrick1","gold1","chartreuse1","dodgerblue1") # c("firebrick2","gold2","chartreuse2","deepskyblue2","slateblue2")
COLS <- rep("firebrick1",length(COMPILE$VAR)) # COLS.list # rep( COLS.list[1], nrow(COMPILE$VAR) )
COLS[ grep("SJC",names(COMPILE$VAR)) ] <- "gold4"
COLS[ grep("TJC",names(COMPILE$VAR)) ] <- "dodgerblue4"
COLS[ grep("CRP",names(COMPILE$VAR)) ] <- "chartreuse4"
COLS[ grep("rSJC",names(COMPILE$VAR)) ] <- "gold1"
COLS[ grep("rTJC",names(COMPILE$VAR)) ] <- "dodgerblue1"
COLS[ grep("lCRP",names(COMPILE$VAR)) ] <- "chartreuse1"
LTYS <- rep( 2, length(COMPILE$VAR) )
LTYS[ grep("DEL",rownames(COMPILE$VAR)) ] <- 1
# PCHS <- rep( "o", length(COMPILE$VAR) )
# PCHS[ grep("MNa",names(COMPILE$VAR)) ] <- "A"
# PCHS[ grep("Bwk",names(COMPILE$VAR)) ] <- "B"
# PCHS[ grep("PRC",names(COMPILE$VAR)) ] <- "C"
# PCHS[ grep("MNcd",names(COMPILE$VAR)) ] <- "D"
# PCHS[ grep("VARdr",names(COMPILE$VAR)) ] <- "E"
# PCHS[ grep("VARwk",names(COMPILE$VAR)) ] <- "F"
# PCHS[ grep("POST",names(COMPILE$VAR)) ] <- "G"
PCHS <- rep( "o", length(COMPILE$VAR) )
PCHS[ grep("WAG4",names(COMPILE$VAR)) ] <- "A"
PCHS[ grep("WAG12",names(COMPILE$VAR)) ] <- "B"
PCHS[ grep("WAG20",names(COMPILE$VAR)) ] <- "C"
PCHS[ grep("WAG28",names(COMPILE$VAR)) ] <- "D"
PCHS[ grep("FL",names(COMPILE$VAR)) ] <- "E"



# load( "/Users/kstandis/Data/Burn/Results/20150423_GCTA_Perm_Test.Rdata" )
# Num_Perms <- nrow( COMPILE$MOD[[1]] ) - 1

## Plot LRT Stat Distribution vs Actual Results
P.perm.comp <- numeric( length(COMPILE$MOD) )
names(P.perm.comp) <- names(COMPILE$MOD)
derive <- "" ; count <- 0
for ( p in 1:length(COMPILE$MOD) ) {
	pheno <- names(COMPILE$MOD)[p]
	derive.old <- derive
	derive <- paste(strsplit(names(COMPILE$MOD)[p],"_")[[1]][1:2],collapse="_")
	if ( derive!=derive.old ) {
		if ( p > 1 ) { dev.off() }
		jpeg( paste(PathToData,"/Perm_LRT_Distr_",derive,".jpeg",sep=""), width=1600,height=1600, pointsize=30 )
		par(mfrow=c(3,3))
		count <- 1
	}
	# Plot LRT Distribution
	XLIM <- range( COMPILE$MOD[[p]][,"LRT"] )
	BRKS <- seq( 0,XLIM[2]+.5,.5 )
	hist( COMPILE$MOD[[p]][-(Num_Perms+1),"LRT"], xlim=XLIM,main=paste("LRT Stat:",pheno),xlab="LRT Stat",breaks=BRKS,col=COLS[p] )
	arrows( COMPILE$MOD[[p]]["True","LRT"],.2*Num_Perms,COMPILE$MOD[[p]]["True","LRT"],0, col=gsub("1","4",COLS[p]),lwd=3 )
	# Print P-Values
	P.perm <- (1+length(which(COMPILE$MOD[[p]][1:Num_Perms,"LRT"]>COMPILE$MOD[[p]]["True","LRT"]))) / (1+Num_Perms)
	P.perm.comp[p] <- P.perm
	P.dat <- COMPILE$MOD[[p]]["True","Pval"]
	text( quantile(XLIM,.2),.8*Num_Perms, pos=4,label=paste("P.perm =",formatC(P.perm,digits=2,format="e")) )
	text( quantile(XLIM,.2),.75*Num_Perms, pos=4,label=paste("P.dat =",formatC(P.dat,digits=2,format="e")) )
	count <- count + 1
	# if ( count>9 ) { dev.off() }
} ; dev.off()

## Plot Permuted vs Actual P-Values
P.dat.comp <- unlist(lapply( COMPILE$MOD, function(x) x["True","Pval"] ))
PCHS.leg <- data.frame( LABS=unlist(lapply(strsplit(names(P.dat.comp),"_"),function(x) paste(x[1],x[2],sep="_"))), PCHS )
PCHS.leg <- PCHS.leg[ which(!duplicated(PCHS.leg$PCHS)), ]
COLS.leg <- data.frame( LABS=unlist(lapply(strsplit(names(P.dat.comp),"_"),function(x) x[3])), COLS )
COLS.leg <- COLS.leg[ which(!duplicated(COLS.leg$COLS)), ]
LIM <- c(0, max(-log10( c(P.dat.comp,P.perm.comp) )) )
jpeg( paste(PathToData,"/Perm_Pvals.jpeg",sep=""), width=1600,height=1600, pointsize=30 )
plot( 0,0, type="n", xlim=LIM,ylim=LIM, xlab="GCTA: -log10(p)",ylab="Permuted: -log10(p)", main="Permuted vs GCTA P-Values" )
abline( h=seq(0,LIM[2]+1,1),lty=2,col="grey50") ; abline( v=seq(0,LIM[2]+1,1),lty=2,col="grey50")
abline( 0,1, lty=1,col="black",lwd=2 )
abline( h=-log10(.05), lty=2,col="magenta2",lwd=3 )
abline( v=-log10(.05), lty=2,col="magenta2",lwd=3 )
abline( h=-log10(1/(1+Num_Perms)), lty=2,col="chocolate2",lwd=3 )
points( -log10(P.dat.comp), -log10(P.perm.comp), col=COLS,pch=PCHS,cex=2 )
legend( "topleft", legend=PCHS.leg$LABS[order(PCHS.leg$PCHS)], pch=as.character(PCHS.leg$PCHS[order(PCHS.leg$PCHS)]), col="black", cex=1.5 )
legend( "bottomright", legend=COLS.leg$LABS[order(COLS.leg$COLS)], col=as.character(COLS.leg$COLS[order(COLS.leg$COLS)]), pch=20, cex=1.5, ncol=2 )
dev.off()

# FRAME <- data.frame( GCTA=P.dat.comp, PERM=P.perm.comp, PCHS, COLS )
# # Plot Confidence Intervals of Permuted and Actual Data
# XLIM <- c( 0,length(COMPILE$MOD)+1 )
# YLIM <- c( -.5,1.5 )
# # COLS <- c("firebrick1","gold1","chartreuse1","dodgerblue1","chartreuse1","dodgerblue1") 
# plot( 0,0,type="n", ylim=YLIM,xlim=XLIM, main="Heritability Estimates for True/Permuted Data", ylab="Heritability Estimate (%)", xlab="Phenotype",xaxt="n" )
# axis( 1, at=1:length(COMPILE$MOD), label=names(COMPILE$MOD),las=2 )
# abline( h=seq(-2,2,.2),lty=2,col="grey50" )
# abline( h=c(0,1),lty=1,col="black" )
# for ( p in 1:length(COMPILE$MOD) ) {
# 	arrows( p, COMPILE$VAR[[p]][,"VgVp"]+COMPILE$SE[[p]][,"VgVp"], p, COMPILE$VAR[[p]][,"VgVp"]-COMPILE$SE[[p]][,"VgVp"], code=3,angle=90, col=COLS[p],lwd=1 )
# }
# for ( p in 1:length(COMPILE$MOD) ) {
# 	arrows( p, COMPILE$VAR[[p]][Num_Perms+1,"VgVp"]+COMPILE$SE[[p]][Num_Perms+1,"VgVp"], p, COMPILE$VAR[[p]][Num_Perms+1,"VgVp"]-COMPILE$SE[[p]][Num_Perms+1,"VgVp"], code=3,angle=90, col=gsub("1","4",COLS[p]),lwd=3 )
# 	points( p, COMPILE$VAR[[p]][Num_Perms+1,"VgVp"], col=gsub("1","4",COLS[p]),pch=20 )
# }



















