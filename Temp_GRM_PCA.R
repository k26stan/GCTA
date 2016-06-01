## Compare GRMs from SNPs at MAF 1 & 5 ##
library(gplots)

## Assign Function to Load Binary Genetic Relationship Files
ReadGRMBin <- function(prefix, AllN=F, size=4) {
	sum_i <- function(i) {
		return(sum(1:i))
	}
	BinFileName <- paste(prefix,".grm.bin",sep="")
	NFileName <- paste(prefix,".grm.N.bin",sep="")
	IDFileName <- paste(prefix,".grm.id",sep="")
	id <- read.table(IDFileName)
	n <- dim(id)[1]
	BinFile <- file(BinFileName, "rb")
	grm <- readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
	NFile <- file(NFileName, "rb")
	if(AllN==T){
		print( "AllN == T" )
		N <- readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
	}else{
		print( "AllN != T" )
		N <- readBin(NFile, n=1, what=numeric(0), size=size)
	}
	i <- sapply(1:n, sum_i)
	closeAllConnections()
	return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}


## Loop Through and Load Genetic Relationship Estimates
PathToMAF.1 <- "/projects/janssen/Heritability/20160125_GCTA_Derived_MAF1_SNP_Manu_PhenoCovs_Derived/1_GRM/1-GRM_FULL"
PathToMAF.5 <- "/projects/janssen/Heritability/20160125_GCTA_Derived_MAF5_SNP_Manu_PhenoCovs_Derived/1_GRM/1-GRM_FULL"
PathToMAF.EUR <- "/projects/janssen/Heritability/20160126_GCTA_Derived_MAF1_SNP_EUR_Manu_PhenoCovs_Derived/1_GRM/1-GRM_FULL"
AllN <- T
size <- 4
DAT.1 <- ReadGRMBin( PathToMAF.1, AllN, size )
DAT.5 <- ReadGRMBin( PathToMAF.5, AllN, size )
DAT.E <- ReadGRMBin( PathToMAF.EUR, AllN, size )

DAT.LIST <- list( DAT.1, DAT.5, DAT.E )
names(DAT.LIST) <- c("MAF.1","MAF.5","EUR")

#################################################
## REFORMAT & PLOT DATA #########################
#################################################
DAT_ARR <- list()
for ( d in 1:length(DAT.LIST) ) {
	DAT.name <- names(DAT.LIST)[d]
	DAT <- DAT.LIST[[d]]
	file_name <- names(DAT.LIST)[d]

	## Set Diagonal Matrix and then Fill
	TEMP_ARR <- diag(DAT$diag) # * DAT$diag
	TEMP_ARR[which(upper.tri(TEMP_ARR))] <- DAT$off
	TEMP_ARR[which(lower.tri(TEMP_ARR))] <- t( TEMP_ARR )[which(lower.tri(TEMP_ARR))]
	rownames(TEMP_ARR) <- colnames(TEMP_ARR) <- DAT$id[,1]
	DAT_ARR[[DAT.name]] <- TEMP_ARR
}


CUT <- function( cutoff ) {
	YUP <- list()
	for ( d in 1:length(DAT.LIST) ) {
		DAT.name <- names(DAT.LIST)[d]
		N <- nrow(DAT_ARR[[DAT.name]])
		YUP[[DAT.name]] <- unlist(lapply( 1:N, function(x) any(DAT_ARR[[DAT.name]][-x,x] > cutoff) ))
	}
	HOW_MANY <- lapply( YUP, function(x)length(which(x==T)) )

	return(HOW_MANY)
}
CUT( .1 )

#################################################
## IDENTIFY BATCH EFFECTS #######################
#################################################

## Load Group Designations
Group_Names <- as.character( read.table( "/projects/janssen/scripts/GROUP_NAMES_ALL.txt" )[,1] )
Group_Names <- Group_Names[ order( as.numeric(sapply(strsplit(Group_Names,"_"),"[",2)) )]
PathToGrps <- "/projects/janssen/scripts/SAMP_GRPS/"
GRPS <- list()
for ( group in Group_Names ) {
	GRPS[[group]] <- as.character( read.table(paste(PathToGrps,group,sep=""))[,1] )
}

## Load Phenotype Table
PathToPheno <- "/projects/janssen/ASSOCIATION/PH-PHENOTYPES/20150520_Full_Table.txt"
FT <- read.table( PathToPheno, sep="\t", header=T )
ANC_COLS <- c("firebrick2","chocolate2","gold1","springgreen2","steelblue2","slateblue3")
jpeg( paste("/home/kstandis/20160125_GRM.Admix.jpeg",sep=""),height=1600,width=2400,pointsize=24 )
par(mfrow=c(4,5))
for ( group in Group_Names ) {
	samps <- GRPS[[group]]
	TEMP_ARR <- FT[ which(as.character(FT$ID) %in% samps), c("europe","africa","amerind","eastAsian","oceania","centralAsia") ]
	rownames(TEMP_ARR) <- FT[ which(as.character(FT$ID) %in% samps),"ID" ]
	barplot( t(TEMP_ARR), col=ANC_COLS, main=group )
}
dev.off()

## FCT: Plot a particular GRM
PLOT <- function( array, tag, FT ) {
	## GRM Plot Parameters
	COLORS <- c("firebrick2","sienna1","yellow1","black","chartreuse1","deepskyblue1","slateblue2")
	N_COLS <- 100
	COLS <- colorRampPalette(COLORS)(N_COLS)
	BRKS <- seq( -1,1,length.out=N_COLS+1)
	JNJ_ANC_COLS <- c("firebrick2","mediumpurple2","gold1","chartreuse1")[factor(FT$JNJ_ANC)] # ,"deepskyblue1")
	ETHN_COLS <- c("chocolate2","deepskyblue1")[factor(FT$ETHN)]
	# MG <- merge( FT[,c("ID","JNJ_ANC","ETHN")], array, by.x="ID",by.y="row.names")
	MG <- merge( data.frame(ID=FT$ID,JNJ_ANC_COLS,ETHN_COLS), array, by.x="ID",by.y="row.names")
	MG <- MG[,c("ID","JNJ_ANC_COLS","ETHN_COLS",as.character(MG$ID))]
	## Plot GRMs
	jpeg( paste("/home/kstandis/20160125_GRM.",tag,".jpeg",sep=""),height=1600,width=1600,pointsize=24 )
	# heatmap.2(data.matrix(MG[,-c(1:3)]), main=paste("Genetic Relationship Matrix -",tag), scale="none", symm=F, dendrogram="both", Rowv=T, Colv=T, col=COLS, breaks=BRKS, sepwidth=c(0,0), trace="none",ColSideColors=JNJ_ANC_COLS[factor(MG$JNJ_ANC)],RowSideColors=ETHN_COLS[factor(MG$ETHN)] )
	heatmap.2(data.matrix(MG[,-c(1:3)]), main=paste("Genetic Relationship Matrix -",tag), scale="none", symm=F, dendrogram="both", Rowv=T, Colv=T, col=COLS, breaks=BRKS, sepwidth=c(0,0), trace="none",ColSideColors=as.character(MG$JNJ_ANC_COLS),RowSideColors=as.character(MG$ETHN_COLS) )
	dev.off()
}

## Plot GRM for each Group at various MAFs
for ( d in 1:length(DAT.LIST) ) {
	DAT.name <- names(DAT.LIST)[d]
	for ( group in Group_Names ) {
		tag <- paste( group,DAT.name,sep="-" )
		samps.1 <- GRPS[[group]]
		array <- DAT_ARR[[DAT.name]]
		colnames(array) <- rownames(array) <- sapply(strsplit(colnames(array),"-"),"[",1)
		samps <- intersect( samps.1, rownames(array) )
		if ( length(samps)>1 ) {
			array <- array[samps,samps]
			PLOT( array, tag, FT )	
		}
	}
}
	

## Plot GRM for Europeans
EUR.samps.old <- as.character( read.table("/projects/janssen/VCFs/PLINK/EUR_SAMPS.txt",header=F)[,1] )
EUR.samps.prc <- as.character( FT$ID_2[which(FT$europe>=.90)] )
EUR.samps.self <- as.character( FT$ID_2[which(is.na(FT$europe) & FT$RACE=="WHITE" & FT$ETHN=="NOT HISPANIC OR LATINO" )] )
# EUR.samps <- as.character( FT$ID[which( FT$RACE=="WHITE" & FT$ETHN=="NOT HISPANIC OR LATINO" )] )
# EUR.samps <- as.character( FT$ID[which(FT$europe>=.95 | (FT$RACE=="WHITE" & FT$ETHN=="NOT HISPANIC OR LATINO" ) )])
TEMP <- FT[, c("COUN","europe","RACE","ETHN","JNJ_ANC","africa") ] ; rownames(TEMP) <- FT[,"ID_2"]
EUR.samps <- union( EUR.samps.prc, EUR.samps.self)
setdiff( EUR.samps.old, EUR.samps )
FT[ which( as.character(FT$ID_2) %in% EUR.samps), c("ID","COUN","europe","RACE","ETHN","JNJ_ANC","africa") ]
FT[ which( !(as.character(FT$ID_2) %in% EUR.samps)), c("ID","COUN","europe","RACE","ETHN","JNJ_ANC","africa") ]
write.table( EUR.samps, "/projects/janssen/Heritability/EUR.samps.txt", col.names=F,row.names=F,quote=F)

## Plot Admixture Estimate Heatmaps
ADMIX_COLS <- c("firebrick2","sienna1","yellow1","chartreuse1","deepskyblue1","slateblue2","black")
N_COLS <- 100
COLS <- colorRampPalette(rev(ADMIX_COLS))(N_COLS)
BRKS <- seq(0,1,length.out=N_COLS+1)
JNJ_ANC_COLS <- c("firebrick2","mediumpurple2","gold1","chartreuse1")[factor(FT$JNJ_ANC)] # ,"deepskyblue1")
ETHN_COLS <- c("chocolate2","deepskyblue1")[factor(FT$ETHN)]
TEMP_ARR <- data.matrix(FT[,c("europe","africa","amerind","eastAsian","oceania","centralAsia")])
rownames(TEMP_ARR) <- as.character(FT$ID)
ORDER.1 <- order( TEMP_ARR[,"eastAsian"] )
ORDER.2 <- order( TEMP_ARR[ORDER.1,"amerind"] )
ORDER.3 <- order( TEMP_ARR[ORDER.1[ORDER.2],"europe"] )
ORDER <- ORDER.1[ORDER.2[ORDER.3]]
jpeg( paste("/home/kstandis/20160125_GRM.Admix.Heat.jpeg",sep=""),height=1600,width=1600,pointsize=24 )
heatmap.2(TEMP_ARR[ORDER,], main=paste("Genetic Relationship Matrix -",tag), scale="none", symm=F, dendrogram="none", Rowv=F, Colv=F, col=COLS, breaks=BRKS, sepwidth=c(0,0), trace="none",RowSideColors=JNJ_ANC_COLS[ORDER] )
dev.off()


PC <- prcomp( na.omit(TEMP_ARR), scale=T )
# PC <- princomp( na.omit(TEMP_ARR), scale=F )
jpeg( paste("/home/kstandis/20160125_GRM.Admix.PC1.jpeg",sep=""),height=1600,width=1600,pointsize=24 )
pairs( PC$x, col=JNJ_ANC_COLS,pch=as.numeric(factor(FT$ETHN)) )
dev.off()
jpeg( paste("/home/kstandis/20160125_GRM.Admix.PC2.jpeg",sep=""),height=800,width=1600,pointsize=24 )
barplot( t(t(PC$rotation)*PC$sdev), beside=T, legend=T )
dev.off()





TEMP_ARR.2 <- data.frame( TEMP_ARR, TEST=TEMP_ARR[,6]+rnorm(nrow(TEMP_ARR),0,.02) )
PC <- princomp( na.omit(TEMP_ARR.2), scale=F )
PC <- prcomp( na.omit(TEMP_ARR.2), scale=T )
