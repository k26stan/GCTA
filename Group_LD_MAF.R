## Compile Variants by MAF/LD for GCTA ##
## September 4, 2015 ##
## Kristopher Standish ##

## After calculating LD for each chromosome, we need to combine them and
 # then group them by LD quantile and MAF.
## Output is list of SNPs for each LD/MAF group
 # PathToSave = PathToLD/*.GRP.txt

###################################################
## PARSE COMMAND LINE #############################
###################################################

LINE <- commandArgs(trailingOnly = TRUE)
# LINE <- "/projects/janssen/HLA"
# LINE <- "/Users/kstandis/Downloads/20150904_GCTA_Test/ART3001"
PathToLD <- LINE[1]
Num_LD_Grps <- LINE[2]

print( "Running: Group_LD_MAF.R" )

###################################################
## LOAD DATA ######################################
###################################################

## Load LD Files
print( "Loading LD Files" )
File_List <- grep( "ld$", list.files(PathToLD), value=T )
LD <- list()
for ( file in File_List ) {
	tag <- gsub( "0-LD_","", gsub(".mrsq.ld","",file, fixed=T) )
	LD[[tag]] <- read.table( paste(PathToLD,file,sep="/"),header=T, stringsAsFactors=F )
	print(paste( tag,"Done" ))
}

## Combine Chromosomes into Single Group
print( "Merging LD Files" )
LD.full <- Reduce( rbind, LD )
LD.full <- data.frame( LD.full, LD_SCORE=LD.full$mean_rsq*LD.full$snp_num )

###################################################
## GROUP BY MAF & LD ##############################
###################################################
print( "Grouping Variants by LD/MAF" )

## Create Keys for MAF & LD Groups
# RNG.MAF <- c( 0, .001, .01, .05, .1, .25, 1 )
RNG.MAF <- c( 0, .01, .05, .1, .25, 1 )
KEY.MAF <- data.frame( CUT=RNG.MAF, GRP=paste("M",RNG.MAF,sep="") )
# QNT.LD <- c( 0, .25, .5, .75, 1 )
QNT.LD <- seq( 0, 1, length.out=Num_LD_Grps+1 )
KEY.LD <- data.frame( CUT=quantile(LD.full$LD_SCORE,QNT.LD), GRP=paste("LDq",0:4,sep="") )

## Group by Level of MAF/LD
GRP.MAF <- cut( LD.full$freq, KEY.MAF$CUT, KEY.MAF$GRP[-1], include.lowest=T )
GRP.LD <- cut( LD.full$LD_SCORE, KEY.LD$CUT, KEY.LD$GRP[-1], include.lowest=T )

## Designate Actual Group based on BOTH
GRP <- paste( GRP.MAF, GRP.LD, sep="_" )
GRP.uniq <- sort(unique(GRP))

###################################################
## SAVE SNP LISTS by GROUP ########################
###################################################
print( "Saving SNP Files" )

## Save Files w/ SNP Lists
for ( grp in GRP.uniq ) {
	PathToSave <- paste(PathToLD,"/",grp,".GRP.txt",sep="")
	SNPs <- LD.full$SNP[ GRP==grp ]
	write.table( SNPs, file=PathToSave, col.names=F,row.names=F,quote=F )
}

###################################################
## JUST FOR KICKS, PLOT LD BY CHROMOSOME ##########
###################################################
print( "Making LD Plot" )

COLS <- c( "dodgerblue4", "tomato1" )
jpeg( paste(PathToLD,"/LD_Plot.jpeg",sep=""), height=2000,width=3000,pointsize=30 )
par(mfrow=c(4,6))
for ( tag in names(LD) ) {
	DAT <- LD[[tag]]
	plot( mean_rsq*snp_num ~ bp, data=DAT, col=COLS[1],pch=20,cex=.5, main=paste("LD:",tag),xlab="Position",ylab="LD_Score" )
	points( mean_lds ~ bp, data=DAT, col=COLS[2],type="l",lwd=2 )
}
dev.off()

###################################################
## END OF DOC #####################################
###################################################
