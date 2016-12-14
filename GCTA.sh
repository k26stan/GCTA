## Run GCTA on Janssen Psych Cohorts ##
## Based on Nathan's scripts ##
## January 26, 2015 ##
## Kristopher Standish ##

##########################################################################
## 1 ## Set up Paths #####################################################
##########################################################################
 # Use Bash
 # Take in arguments, set up directories/paths for files/tools
echo \### 1 - `date` \###
echo \### Define Set Variables and Paths \###

###########################################################
## Manually Input Parameters ##

## Names/Paths for Output
DATE=$1
HOME_DIR=$2

## Parameters and Files
VAR_FILE=$3
VAR_DIR=$4
PHENO_DIR=$5
PHENO_FILE=$6
PHENO_NAME_LIST_FILE=$7 # Which Phenotypes and Covariates are you using?
PHENO_NAME_LIST_DIR=$8 # Directory for above file
COV_FILE=$9 # Path to Covariate File or "F"
LD_MAF_GRPS=${10} # LD_MAF
GRM_CUTOFF=${11} # Decimal
PC_COUNT=${12} # How many PCs to Include as Covariates?
START_STEP=${13} # Which Step do you want to start on?

###########################################################
## Constant Paths ##

## Public Tools
GCTA=/projects/janssen/Tools/gcta/gcta64
PLINK=/projects/janssen/Tools/plink_linux_x86_64/plink

## Custom Scripts
GROUP_LD_MAF_R=/projects/janssen/Psych/Scripts/GCTA/Group_LD_MAF.R
PULL_COVS=/projects/janssen/Psych/Scripts/GCTA/Pull_Cov_Cols.R
PULL_PHENO=/projects/janssen/Psych/Scripts/GCTA/Pull_Pheno_Col.R
PLOT_GRM=/projects/janssen/Psych/Scripts/GCTA/Plot_GRM.R
PLOT_EST=/projects/janssen/Psych/Scripts/GCTA/Plot_Estimates.R
PERMUTE=/projects/janssen/Psych/Scripts/GCTA/Permute.R
PLOT_PERMUTED=/projects/janssen/Psych/Scripts/GCTA/Plot_Permuted.R

## Set Specific Paths
VAR_PATH=${VAR_DIR}/${VAR_FILE}
COV_PATH=${PHENO_DIR}/${COV_FILE}
PHENO_PATH=${PHENO_DIR}/${PHENO_FILE}
PHENO_NAME_LIST_PATH=${PHENO_NAME_LIST_DIR}/${PHENO_NAME_LIST_FILE}

## Make new folder for Today's adventures
OUT_DIR=${HOME_DIR}/${DATE}_${PHENO_NAME_LIST_FILE%%.txt}
mkdir ${OUT_DIR}
cd ${OUT_DIR}

###########################################################
## Pull some Info out of Parameters ## 

## Determine if I'm using Covariates
if [ -e ${COV_PATH} ]
then
USE_COVARS=TRUE
else
USE_COVARS=FALSE
fi

## Determine if GRM Cutoff is Used
if (( $(echo "$GRM_CUTOFF > 0" | bc -l) )) ; then
USE_GRM=${OUT_DIR}/1_GRM/1-GRM_CUT
else
USE_GRM=${OUT_DIR}/1_GRM/1-GRM_FULL
fi
echo $USE_GRM

## Pull out number of LD & MAF Groups
LD_GRPS=`echo $LD_MAF_GRPS | cut -d_ -f1`
MAF_GRPS=`echo $LD_MAF_GRPS | cut -d_ -f2`

## Specify a File to which to Write Updates
UPDATE_FILE=${OUT_DIR}/Update.txt

## Done
if [ "$START_STEP" -le 1 ]; then
echo `date` "1 - Define Set Variables and Paths - DONE" > ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 2 ## Calculate Linkage Disequilibrium #################################
##########################################################################
if [ "$START_STEP" -le 2 ]; then
echo \### 2 - `date` \###
echo \### Calculate LD \###
echo `date` "2 - Calculate LD" >> ${UPDATE_FILE}

##########################################################
## Prep for LD Calculation ##

## FCT: Calculate minimum MAF
calc() {
awk "BEGIN { print "$*" }"
}

## FCT: Calculate LD by Chromosome
LD_By_Chrom() {
chr_list=$1
for chr in $chr_list; do
 # Pull out Chromosome for LD Analysis
${PLINK} \
--bfile ${VAR_PATH} \
--make-bed \
--maf ${MAF} \
--chr ${chr} \
--out TEMP_CHR_BED.${chr}
 # Calculate LD for Chromosome
${GCTA} \
--bfile TEMP_CHR_BED.${chr} \
--thread-num 1 \
--autosome \
--chr ${chr} \
--ld-score-region 200 \
--out ${OUT_DIR}/0_LD/0-LD_CHR${chr}
 # Remove Temporary BED file
rm TEMP_CHR_BED.${chr}*
done
}

## Calculate Allele Frequency for Removal of Single/Doubletons
N_PATS=`wc -l ${VAR_PATH}.fam | cut -f1 -d " "`
MAF=`calc 1/$N_PATS`
echo $MAF

## Create Path for LD Files
mkdir ${OUT_DIR}/0_LD

## If necessary, Calculate LD Files
if [ "$LD_GRPS" -gt 1 -o "$MAF_GRPS" -gt 1 ] ; then
LD_By_Chrom "`echo {1..22..2}`" &
LD_By_Chrom "`echo {2..22..2}`" &
wait
else
echo "LD Structure Not Calculated"
fi

## Done
echo `date` "2 - Calculate LD - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 3 ## Create LD SNP Groups #############################################
##########################################################################
if [ "$START_STEP" -le 3 ]; then
echo \### 3 - `date` \###
echo \### Create LD Groups \###
echo `date` "3 - Create LD Groups" >> ${UPDATE_FILE}

## If necessary, Determine LD/MAF SNP Groups
if [ "$LD_GRPS" -gt 1 -o "$MAF_GRPS" -gt 1 ] ; then
Rscript ${GROUP_LD_MAF_R} ${OUT_DIR}/0_LD 3
else
echo "LD Structure Not Calculated"
fi

## Done
echo `date` "3 - Create LD Groups - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 4 ## Calculate Genetic Relationship Matrix ############################
##########################################################################
if [ "$START_STEP" -le 4 ]; then
echo \### 4 - `date` \###
echo \### Calculate GRM \###
echo `date` "4 - Calculate GRM" >> ${UPDATE_FILE}

## Create Path for GRM Files
mkdir ${OUT_DIR}/1_GRM

##########################################################
## Calculate GRM(s)
if [ "$LD_GRPS" -gt 1 -o "$MAF_GRPS" -gt 1 ] ; then

 # Using Multiple variants sets from LD/MAF Stratification
for snp_file in `ls ${OUT_DIR}/1_GRM/*GRP.txt`; do
file_name_only=`echo $snp_file | xargs -n1 basename`
${GCTA} \
--bfile ${VAR_PATH} \
--thread-num 1 \
--extract ${snp_file} \
--autosome \
--make-grm \
--out ${OUT_DIR}/1_GRM/1-GRM_FULL.${file_name_only%%.GRP.txt}
done
 # Create List of GRM Files (for --mgrm command)
ls ${OUT_DIR}/1_GRM/* | grep "grm.id" | sed 's/.grm.id//g' > ${OUT_DIR}/GRM_List.txt
 # Put individual GRMs together into 1 (for PCA)
if [ $PC_COUNT -ne 0 ] ; then
${GCTA} \
--mgrm ${OUT_DIR}/GRM_List.txt \
--thread-num 1 \
--make-grm \
--out ${OUT_DIR}/1_GRM/1-GRM_FULL
fi

else
 # Using all Variants at once
${GCTA} \
--bfile ${VAR_PATH} \
--thread-num 1 \
--maf .01 \
--autosome \
--make-grm \
--out ${OUT_DIR}/1_GRM/1-GRM_FULL

fi # Close If LD/MAF Group > 1

## Done
echo `date` "4 - Calculate GRM - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 5 ## Deal w/ Population Structure #####################################
##########################################################################
if [ "$START_STEP" -le 5 ]; then
echo \### 5 - `date` \###
echo \### Population Structure \###
echo `date` "5 - Population Structure" >> ${UPDATE_FILE}

#####################################################
## Do PCA Analysis from GRM
if [ $PC_COUNT -ne 0 ]
then
 ${GCTA} --pca 20 \
 --grm ${OUT_DIR}/1_GRM/1-GRM_FULL \
 --thread-num 1 \
 --out 2-PCA
fi

#########################################################
## Filter GRM relationships above Threshold
if (( $(echo "$GRM_CUTOFF > 0" | bc -l) )) ; then
${GCTA} \
--grm ${OUT_DIR}/1_GRM/1-GRM_FULL \
--thread-num 1 \
--grm-cutoff ${GRM_CUTOFF} \
--make-grm \
--out ${OUT_DIR}/1_GRM/1-GRM_CUT
# USE_GRM=${OUT_DIR}/1_GRM/1-GRM_CUT
# else
# USE_GRM=${OUT_DIR}/1_GRM/1-GRM_FULL
fi

## Done
echo `date` "5 - Population Structure - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi

##########################################################################
## 6 ## Estimate Heritability ############################################
##########################################################################
if [ "$START_STEP" -le 6 ]; then
echo \### 6 - `date` \###
echo \### Estimate Heritability \###
echo `date` "6 - Estimate Heritability" >> ${UPDATE_FILE}

mkdir ${OUT_DIR}/Phenos
mkdir ${OUT_DIR}/3_REML/

## Remove Commands File (if it already exists)
if [ -e ${OUT_DIR}/3_REML/COMMANDS.run ]; then
rm ${OUT_DIR}/3_REML/COMMANDS.run
fi

## Specify Field Separater
IFSo=$IFS
IFS=$'\n' # Makes it so each line is read whole (not separated by tabs)

#####################################################
## Loop through Phenotypes ##########################
# line=`head -1 ${PHENO_NAME_LIST_PATH}`
# line=`head -10 ${PHENO_NAME_LIST_PATH} | tail -1`
for line in `cat ${PHENO_NAME_LIST_PATH}` ; do
echo Next Line: $line
# for line in `head -10 ${PHENO_NAME_LIST_PATH}` ; do

#####################################################
## Pull out Phenotype to Individual File
 # Determine which Phenotype to Use
pheno=`echo ${line} | awk '{print $1}'`
 # Create New Pehnotype File
NEW_PHENO_PATH=${OUT_DIR}/Phenos/${pheno}_FULL.txt
Rscript ${PULL_PHENO} ${PHENO_PATH} ${pheno} ${NEW_PHENO_PATH}

#####################################################
## Specify single or multiple GRMs
if [ "$LD_GRPS" -gt 1 -o "$MAF_GRPS" -gt 1 ] ; then
specify_grm="--mgrm ${OUT_DIR}/GRM_List.txt"
else
specify_grm="--grm ${USE_GRM}"
fi

#####################################################
## Set up Covariate List
 # If Covariates are being Used
if [[ $USE_COVARS == TRUE ]] ; then
 # Determine Which Covariates for this Phenotype (from PHENO_NAME_LIST_FILE)
COVS=`echo ${line} | cut -d$'\t' -f2- | sed 's/\t/QQQ/g'`
 # If PC's specified
if [ $PC_COUNT -gt 0 ] ; then
PCS=`seq 1 ${PC_COUNT}`
PCS_COMMAND=`echo "PC"${PCS} | sed 's/ /QQQPC/g'`
COVS_COMMAND=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}QQQPC${PC_COUNT}" | sed 's/QQQ/_/g'`
else
COVS_COMMAND=`echo "${COVS}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}" | sed 's/QQQ/_/g'`
fi # Close (if PCs)
fi # Close (if USE_COVARS)

#####################################################
## Specify Covariates or Not
if [[ $USE_COVARS == TRUE && $COVS_COMMAND != "" ]] ; then
NEW_COV_PATH=${OUT_DIR}/Phenos/Cov_${COVS_FILENAME}_FULL.txt
 # Compile Covariates into Correct Format
Rscript ${PULL_COVS} ${COVS_COMMAND} ${OUT_DIR}/2-PCA.eigenvec ${COV_PATH} ${NEW_COV_PATH}
specify_cov_pheno="--pheno ${NEW_PHENO_PATH} --qcovar ${NEW_COV_PATH}"
else
specify_cov_pheno="--pheno ${NEW_PHENO_PATH}"
fi

#####################################################
## Set Path for output of GCTA Heritability Estimate
PATH_OUT=${OUT_DIR}/3_REML/3-REML_${pheno}
# echo \## Path to Output: $PATH_OUT
# echo \## GRM flag: $specify_grm
# echo \## Cov/Phenos flag: $specify_cov_pheno

#####################################################
## Compile Commands into .run File
 # Create Command
COMMAND="${GCTA} --reml \
${specify_grm} \
${specify_cov_pheno} \
--reml-maxit 700 \
--reml-est-fix \
--reml-pred-rand \
--out ${PATH_OUT}"
 # Append Commands to "COMMANDS" list file
echo $COMMAND >> ${OUT_DIR}/3_REML/COMMANDS.run

done # Close Phenotype Loop
IFS=$' \t\n' # Reset

## Run Command from List
chmod 777 ${OUT_DIR}/3_REML/COMMANDS.run
${OUT_DIR}/3_REML/COMMANDS.run
# for file in `ls 3_REML/*MNa*hsq`; do echo $file ; grep "V(G)/" $file ; grep "Pval" $file ; done

## Done
echo `date` "6 - Estimate Heritability - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 7 ## Plot GRM & Estimates #############################################
##########################################################################
if [ "$START_STEP" -le 7 ]; then
echo \### 7 - `date` \###
echo \### Make Plots \###
echo `date` "7 - Make Plots" >> ${UPDATE_FILE}

## Plot Genetic Relationship Matrices
Rscript ${PLOT_GRM} ${USE_GRM}

## Plot Heritability Estimates
Rscript ${PLOT_EST} ${PHENO_NAME_LIST_PATH} ${OUT_DIR} ${VAR_FILE}

## Done
echo `date` "7 - Make Plots - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 8 ## Permute for P-Values #############################################
##########################################################################
if [ "$START_STEP" -le 8 ]; then
echo \### 8 - `date` \###
echo \### Permute \###
echo `date` "8 - Permute" >> ${UPDATE_FILE}

mkdir ${OUT_DIR}/4-PERM/
## Specify Number of Permutations
N_PERM=1000

IFSo=$IFS
IFS=$'\n' # Makes it so each line is read whole (not separated by tabs)
## Loop through Phenotypes
# for line in `head -20 ${PHENO_NAME_LIST_PATH}`
for line in `cat ${PHENO_NAME_LIST_PATH}`; do
# for line in `cat ${PHENO_NAME_LIST_PATH} | grep MNcd`; do
# for line in `cat ${PHENO_NAME_LIST_PATH} | grep "PRC\|Bwk\|VARdr\|VARwk"`; do
# for line in `cat ${PHENO_NAME_LIST_PATH} | grep DEL | grep MNa`; do

	# Determine which Phenotype to Use
	pheno=`echo ${line} | awk '{print $1}'`

	## Pull out Phenotype to File
	NEW_PHENO_PATH=${OUT_DIR}/Phenos/${pheno}_FULL.txt

	## If Covariates are being Used
	if [[ $USE_COVARS == TRUE ]]; then
	# Determine Which Covariates for this Phenotype
	# (from PHENO_NAME_LIST_FILE)
	COVS=`echo ${line} | cut -d$'\t' -f2- | sed 's/\t/QQQ/g'`
	# If PC's specified
	if [ $PC_COUNT -eq 0 ]; then
	COVS_COMMAND=`echo "${COVS}" | sed 's/QQQ/,/g'`
	COVS_FILENAME=`echo "${COVS}" | sed 's/QQQ/_/g'`
	else
	PCS=`seq 1 ${PC_COUNT}`
	PCS_COMMAND=`echo "PC"${PCS} | sed 's/ /QQQPC/g'`
	COVS_COMMAND=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/,/g'`
	# COVS_FILENAME=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/_/g'`
	COVS_FILENAME=`echo "${COVS}QQQPC${PC_COUNT}" | sed 's/QQQ/_/g'`
	fi # Close (if PCs)

	## Incorporate Country/Site of Study as Binary Covariate (if Included)
	if [[ $COVS == *COUN* ]]; then
	COVS_COMMAND=`echo $COVS_COMMAND | sed 's/COUN/CN_ARG,CN_AUS,CN_COL,CN_HUN,CN_LTU,CN_MEX,CN_MYS,CN_NZL,CN_POL,CN_RUS,CN_UKR/g'`
	fi # Close (if COUN)

	## Compile Covariates into Correct Format
	NEW_COV_PATH=${OUT_DIR}/Phenos/Cov_${COVS_FILENAME}_FULL.txt
	echo ${NEW_COV_PATH}

	fi # Close (if USE_COVARS)

	#######################################
	## Permute Phenotype/Covariate Files ##
	Rscript ${PERMUTE} ${NEW_PHENO_PATH} ${NEW_COV_PATH} ${N_PERM} ${COVS_COMMAND}

	## Loop Through Z Permutations
	for perm in `seq ${N_PERM}`; do

		## Set up Path for GCTA Output
		EST_OUT=${OUT_DIR}/4-PERM/${pheno}_${perm}

		## Run GCTA to get Heritability Estimates
		# If Covariates are Specified
		if [[ $USE_COVARS == TRUE && $COVS != "" ]]; then
		${GCTA} \
		--grm ${USE_GRM} \
		--pheno ${NEW_PHENO_PATH%%txt}${perm}.txt \
		--qcovar ${NEW_COV_PATH%%txt}${perm}.txt \
		--reml \
		--reml-maxit 1000 \
		--reml-est-fix \
		--reml-pred-rand \
		--out ${EST_OUT}
		else # If Covariates are NOT Specified
		${GCTA} \
		--grm ${USE_GRM} \
		--pheno ${NEW_PHENO_PATH%%txt}${perm}.txt \
		--reml \
		--reml-maxit 1000 \
		--reml-est-fix \
		--reml-pred-rand \
		--out ${EST_OUT}
		fi

		rm ${NEW_PHENO_PATH%%txt}${perm}.txt
		rm ${NEW_COV_PATH%%txt}${perm}.txt

	done # Close Permutation Loop

done # Close Phenotype Loop

IFS=$IFSo # Reset


## Done
echo `date` "8 - Permute - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 9 ## Plot Permuted Results ############################################
##########################################################################
if [ "$START_STEP" -le 9 ]; then
echo \### 9 - `date` \###
echo \### Plot Permuted Results \###
echo `date` "9 - Plot Permuted Results" >> ${UPDATE_FILE}

## Plot Permuted Results
Rscript ${PLOT_PERMUTED} ${PHENO_NAME_LIST_PATH} ${OUT_DIR} ${N_PERM}

## Done
echo `date` "9 - Plot Permuted Results - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## END OF DOC ############################################################
##########################################################################












