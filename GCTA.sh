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
PHENO_NAME_LIST_DIR=$8
COV_FILE=$9 # Path to Covariate File or "F"
PC_COUNT=${10} # How many PCs to Include as Covariates?
START_STEP=${11} # Which Step do you want to start on?

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

# ## Get Names of Specific Files
# TEMP=(${VAR_FILE//\// })
# VAR_FILE_NAME=${TEMP[${#TEMP[@]} - 1]} # Get Name of Variant File

## Determine if I'm using Covariates
if [ -e ${COV_PATH} ]
then
USE_COVARS=TRUE
else
USE_COVARS=FALSE
fi

# ## Specify list of Covariates to include (for command and for filename)
# if [[ $USE_COVARS == TRUE ]]
# then

# if [ $PC_COUNT -eq 0 ]
# then
# COVS_COMMAND=`echo "${COVS}" | sed 's/QQQ/,/g'`
# COVS_FILENAME=`echo "${COVS}" | sed 's/QQQ/_/g'`
# else
# PCS=`seq 1 ${PC_COUNT}`
# PCS_COMMAND=`echo "PC"${PCS} | sed 's/ /QQQPC/g'`
# COVS_COMMAND=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/,/g'`
# COVS_FILENAME=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/_/g'`
# fi

# ## Incorporate Country/Site of Study as Binary Covariate (if Included)
# if [[ $COVS == *COUN* ]]
# then
# COVS_COMMAND=`echo $COVS_COMMAND | sed 's/COUN/CN_ARG,CN_AUS,CN_COL,CN_HUN,CN_LTU,CN_MEX,CN_MYS,CN_NZL,CN_POL,CN_RUS,CN_UKR/g'`
# fi

# fi # Close (if USE_COVARS)
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

## Calculate Allele Frequency
 # For Removal of Single/Doubletons
N_PATS=`wc -l ${VAR_PATH}.fam | cut -f1 -d " "`
calc() {
	awk "BEGIN { print "$*" }"
}
MAF=`calc 1/$N_PATS`
echo $MAF

## Create Path for LD Files
mkdir ${OUT_DIR}/0_LD

## FCT: Calculate LD by Chromosome
LD_By_Chrom() {
chr_list=$1
for chr in $chr_list; do
## Pull out Chromosome for LD Analysis
${PLINK} \
--bfile ${VAR_PATH} \
--make-bed \
--maf ${MAF} \
--chr ${chr} \
--out TEMP_CHR_BED.${chr}
## Calculate LD for Chromosome
${GCTA} \
--bfile TEMP_CHR_BED.${chr} \
--thread-num 1 \
--autosome \
--chr ${chr} \
--ld-score-region 200 \
--out ${OUT_DIR}/0_LD/0-LD_CHR${chr}
## Remove Temporary BED file
rm TEMP_CHR_BED.${chr}*
done
}

## Run Function for LD Calculation
LD_By_Chrom "`echo {1..22..2}`" &
LD_By_Chrom "`echo {2..22..2}`" &
wait

## Done
echo `date` "2 - Calculate LD - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 3 ## Calculate Genetic Relationship Matrix ############################
##########################################################################
if [ "$START_STEP" -le 3 ]; then
echo \### 3 - `date` \###
echo \### Calculate GRM \###
echo `date` "3 - Calculate GRM" >> ${UPDATE_FILE}

##########################################################
## Group Variants by LD & MAF
Rscript ${GROUP_LD_MAF_R} ${OUT_DIR}/0_LD 3

## Create Path for GRM Files
mkdir ${OUT_DIR}/1_GRM

##########################################################
## Calculate GRM Using all variants
for snp_file in `ls ${OUT_DIR}/0_LD/*GRP.txt`; do
## Pull out Group of SNPs for GRM Calculation
${PLINK} \
--bfile ${VAR_PATH} \
--make-bed \
--extract ${snp_file} \
--out TEMP_BED
## Calculate GRM on Group of SNPs
file_name_only=`echo $snp_file | xargs -n1 basename`
${GCTA} \
--bfile TEMP_BED \
--thread-num 1 \
--autosome \
--make-grm \
--out ${OUT_DIR}/1_GRM/1-GRM_FULL.${file_name_only%%.GRP.txt}
done

## Remove Temp BED File
rm TEMP_BED*
## Create List of GRM Files (for --mgrm command)
ls ${OUT_DIR}/1_GRM/* | grep "grm.id" | sed 's/.grm.id//g' > ${OUT_DIR}/GRM_List.txt

##########################################################
## Filter relationships below .05
# ${GCTA} \
# --grm 1-GRM_FULL \
# --thread-num 1 \
# --grm-cutoff 0.05 \
# --make-grm \
# --out 1-GRM_FULL.RM5
# for file in `ls ${OUT_DIR}/1_GRM/1-GRM_FULL*`
# do
# extension=`echo ${file} | sed 's/1-GRM_FULL//g'`
# new_file=1-GRM_FULL.RM5${extension}
# cp ${file} ${new_file}
# done

## Done
echo `date` "3 - Calculate GRM - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
# ##########################################################################
# ## 3 ## Calculate Principal Components ###################################
# ##########################################################################
# if [ "$START_STEP" -le 3 ]; then
# echo \### 3 - `date` \###
# echo \### Calculate PCs \###
# echo `date` "3 - Calculate PCs" >> ${UPDATE_FILE}

# #####################################################
# ## Do PCA Analysis for each GRM
# for snp_file in `ls ${OUT_DIR}/0_LD/*GRP.txt`; do
#  ${GCTA} --pca 20 \
#  --grm ${OUT_DIR}/1_GRM/1-GRM_FULL.RM5 \
#  --thread-num 1 \
#  --out 2-PCA_FULL

# ## Done
# echo `date` "3 - Calculate PCs - DONE" >> ${UPDATE_FILE}
# printf "V\nV\nV\nV\nV\nV\nV\nV\n"
# fi
##########################################################################
## 4 ## Pull out Covariates ##############################################
##########################################################################
# if [ "$START_STEP" -le 4 ]; then
# echo \### 4 - `date` \###
# echo \### Collect Covariates \###
# echo `date` "4 - Collect Covariates" >> ${UPDATE_FILE}

# ## If Covariates are being Used
# if [[ $USE_COVARS == TRUE ]]
# then

# COVS=`
# ## Specify list of Covariates to include (for command and for filename)
# if [[ $USE_COVARS == TRUE ]]
# then

# if [ $PC_COUNT -eq 0 ]
# then
# COVS_COMMAND=`echo "${COVS}" | sed 's/QQQ/,/g'`
# COVS_FILENAME=`echo "${COVS}" | sed 's/QQQ/_/g'`
# else
# PCS=`seq 1 ${PC_COUNT}`
# PCS_COMMAND=`echo "PC"${PCS} | sed 's/ /QQQPC/g'`
# COVS_COMMAND=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/,/g'`
# COVS_FILENAME=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/_/g'`
# fi

# ## Incorporate Country/Site of Study as Binary Covariate (if Included)
# if [[ $COVS == *COUN* ]]
# then
# COVS_COMMAND=`echo $COVS_COMMAND | sed 's/COUN/CN_ARG,CN_AUS,CN_COL,CN_HUN,CN_LTU,CN_MEX,CN_MYS,CN_NZL,CN_POL,CN_RUS,CN_UKR/g'`
# fi

# fi # Close (if USE_COVARS)

# #####################################################
# ## Compile Covariates into Correct Format
# NEW_COV_PATH=${OUT_DIR}/${COVS_FILENAME}_FULL.txt
# echo ${NEW_COV_PATH}
# Rscript ${PULL_COVS} ${COVS_COMMAND} ${OUT_DIR}/2-PCA_FULL.eigenvec ${COV_PATH} ${NEW_COV_PATH}

# fi

# ## Done
# echo `date` "4 - Collect Covariates - DONE" >> ${UPDATE_FILE}
# printf "V\nV\nV\nV\nV\nV\nV\nV\n"
# fi
##########################################################################
## 5 ## Estimate Heritability ############################################
##########################################################################
if [ "$START_STEP" -le 5 ]; then
echo \### 5 - `date` \###
echo \### Estimate Heritability \###
echo `date` "5 - Estimate Heritability" >> ${UPDATE_FILE}

mkdir ${OUT_DIR}/Phenos
mkdir ${OUT_DIR}/3_REML/

IFSo=$IFS
IFS=$'\n' # Makes it so each line is read whole (not separated by tabs)
## Loop through Phenotypes
# for line in `head -20 ${PHENO_NAME_LIST_PATH}`
for line in `cat ${PHENO_NAME_LIST_PATH}`
do
 # Determine which Phenotype to Use
pheno=`echo ${line} | awk '{print $1}'`

 ## If Covariates are being Used
if [[ $USE_COVARS == TRUE ]]
then
 # Determine Which Covariates for this Phenotype
   # (from PHENO_NAME_LIST_FILE)
COVS=`echo ${line} | cut -d$'\t' -f2- | sed 's/\t/QQQ/g'`
 # If PC's specified
if [ $PC_COUNT -eq 0 ]
then
COVS_COMMAND=`echo "${COVS}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}" | sed 's/QQQ/_/g'`
else
PCS=`seq 1 ${PC_COUNT}`
PCS_COMMAND=`echo "PC"${PCS} | sed 's/ /QQQPC/g'`
COVS_COMMAND=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/_/g'`
fi # Close (if PCs)

## Incorporate Country/Site of Study as Binary Covariate (if Included)
if [[ $COVS == *COUN* ]]
then
COVS_COMMAND=`echo $COVS_COMMAND | sed 's/COUN/CN_ARG,CN_AUS,CN_COL,CN_HUN,CN_LTU,CN_MEX,CN_MYS,CN_NZL,CN_POL,CN_RUS,CN_UKR/g'`
fi # Close (if COUN)

fi # Close (if USE_COVARS)

## Pull out Phenotype to File
NEW_PHENO_PATH=${OUT_DIR}/Phenos/${pheno}_FULL.txt
Rscript ${PULL_PHENO} ${PHENO_PATH} ${pheno} ${NEW_PHENO_PATH}

EST_OUT=${OUT_DIR}/3_REML/3-REML_${pheno}

## Run GCTA to get Heritability Estimates
 # If Covariates are Specified
if [[ $USE_COVARS == TRUE && $COVS != "" ]]
then
## Compile Covariates into Correct Format
NEW_COV_PATH=${OUT_DIR}/Phenos/${COVS_FILENAME}_FULL.txt
echo ${NEW_COV_PATH}
Rscript ${PULL_COVS} ${COVS_COMMAND} ${OUT_DIR}/2-PCA_FULL.eigenvec ${COV_PATH} ${NEW_COV_PATH}
## Run GCTA w/ Covariates
 # gcta64 --reml --mgrm multi_GRMs.txt --pheno phen.txt --out test
${GCTA} \
--reml \
--mgrm ${OUT_DIR}/GRM_List.txt \
--pheno ${NEW_PHENO_PATH} \
--qcovar ${NEW_COV_PATH} \
--reml-maxit 700 \
--reml-est-fix \
--reml-pred-rand \
--out ${EST_OUT}
else
${GCTA} \
--reml \
--mgrm ${OUT_DIR}/GRM_List.txt \
--pheno ${NEW_PHENO_PATH} \
--reml-maxit 700 \
--reml-est-fix \
--reml-pred-rand \
--out ${EST_OUT}
fi

done # Close Phenotype Loop

IFS=$IFSo # Reset

## Done
echo `date` "5 - Estimate Heritability - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 6 ## Plot GRM & Estimates #############################################
##########################################################################
if [ "$START_STEP" -le 6 ]; then
echo \### 6 - `date` \###
echo \### Make Plots \###
echo `date` "6 - Make Plots" >> ${UPDATE_FILE}

## Plot Genetic Relationship Matrices
Rscript ${PLOT_GRM} 1-GRM_FULL

## Plot Heritability Estimates
Rscript ${PLOT_EST} ${PHENO_NAME_LIST_PATH} ${OUT_DIR} ${VAR_FILE}

## Done
echo `date` "6 - Make Plots - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 7 ## Permute for P-Values #############################################
##########################################################################
if [ "$START_STEP" -le 7 ]; then
echo \### 7 - `date` \###
echo \### Permute \###
echo `date` "7 - Permute" >> ${UPDATE_FILE}

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
COVS_FILENAME=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/_/g'`
fi # Close (if PCs)

## Incorporate Country/Site of Study as Binary Covariate (if Included)
if [[ $COVS == *COUN* ]]; then
COVS_COMMAND=`echo $COVS_COMMAND | sed 's/COUN/CN_ARG,CN_AUS,CN_COL,CN_HUN,CN_LTU,CN_MEX,CN_MYS,CN_NZL,CN_POL,CN_RUS,CN_UKR/g'`
fi # Close (if COUN)

## Compile Covariates into Correct Format
NEW_COV_PATH=${OUT_DIR}/Phenos/${COVS_FILENAME}_FULL.txt
echo ${NEW_COV_PATH}

fi # Close (if USE_COVARS)

#######################################
## Permute Phenotype/Covariate Files ##
Rscript ${PERMUTE} ${NEW_PHENO_PATH} ${NEW_COV_PATH} ${N_PERM}

## Loop Through Z Permutations
for perm in `seq ${N_PERM}`; do

## Set up Path for GCTA Output
EST_OUT=${OUT_DIR}/4-PERM/${pheno}_${perm}

## Run GCTA to get Heritability Estimates
# If Covariates are Specified
if [[ $USE_COVARS == TRUE && $COVS != "" ]]; then
${GCTA} \
--grm 1-GRM_FULL.RM5 \
--pheno ${NEW_PHENO_PATH%%txt}${perm}.txt \
--qcovar ${NEW_COV_PATH%%txt}${perm}.txt \
--reml \
--reml-maxit 1000 \
--reml-est-fix \
--reml-pred-rand \
--out ${EST_OUT}
else # If Covariates are NOT Specified
${GCTA} \
--grm 1-GRM_FULL.RM5 \
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
echo `date` "7 - Permute - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 8 ## Plot Permuted Results ############################################
##########################################################################
if [ "$START_STEP" -le 8 ]; then
echo \### 8 - `date` \###
echo \### Plot Permuted Results \###
echo `date` "8 - Plot Permuted Results" >> ${UPDATE_FILE}

## Plot Permuted Results
Rscript ${PLOT_PERMUTED} ${PHENO_NAME_LIST_PATH} ${OUT_DIR} ${N_PERM}

## Done
echo `date` "8 - Plot Permuted Results - DONE" >> ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## END OF DOC ############################################################
##########################################################################












