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
PHENO_DIR=$4
PHENO_FILE_LIST=$5 # Which Phenotype Files are you using?
COV_FILE=$6 # Path to Covariate File or "F"
COVS=$7 # Which Covariates to Include?
PC_COUNT=$8 # How many PCs to Include as Covariates?
START_STEP=$9 # Which Step do you want to start on?

###########################################################
## Constant Paths ##

## Public Tools
GCTA=/projects/janssen/Tools/gcta/gcta64 # /gpfs/group/schork/nwineing/gcta/gcta64

## Custom Scripts
PULL_COVS=/projects/janssen/Psych/Scripts/GCTA/Pull_Cov_Cols.R

## Specify Directories
VAR_DIR=/projects/janssen/Psych/Data/Genotyped/

## Set Specific Paths
COV_PATH=${PHENO_DIR}/${COV_FILE}

## Make new folder for Today's adventures
OUT_DIR=${HOME_DIR}/${DATE}_${PHENO_FILE_LIST%%.txt}
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

## Specify list of Covariates to include (for command and for filename)
if [[ $USE_COVARS == TRUE ]]
then

if [ $PC_COUNT -eq 0 ]
then
COVS_COMMAND=`echo "${COVS}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}" | sed 's/QQQ/_/g'`
else
PCS=`seq 1 ${PC_COUNT}`
PCS_COMMAND=`echo "PC"${PCS} | sed 's/ /QQQPC/g'`
COVS_COMMAND=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/,/g'`
COVS_FILENAME=`echo "${COVS}QQQ${PCS_COMMAND}" | sed 's/QQQ/_/g'`
fi

## Incorporate Country/Site of Study as Binary Covariate (if Included)
if [[ $COVS == *COUN* ]]
then
COVS_COMMAND=`echo $COVS_COMMAND | sed 's/COUN/CN_ARG,CN_AUS,CN_COL,CN_HUN,CN_LTU,CN_MEX,CN_MYS,CN_NZL,CN_POL,CN_RUS,CN_UKR/g'`
fi

fi # Close (if USE_COVARS)
## Specify a File to which to Write Updates
UPDATE_FILE=${OUT_DIR}/Update.txt

## Done
if [ "$START_STEP" -le 1 ]; then
echo `date` "1 - Define Set Variables and Paths - DONE" > ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 2 ## Calculate Genetic Relationship Matrix ############################
##########################################################################
if [ "$START_STEP" -le 2 ]; then
echo \### 2 - `date` \###
echo \### Calculate GRM \###
echo `date` "2 - Calculate GRM" >> ${UPDATE_FILE}

##########################################################
## Calculate GRM Using all variants
${GCTA} \
--bfile ${VAR_FILE} \
--thread-num 1 \
--autosome \
--make-grm \
--out 1-GRM_FULL

##########################################################
## Filter relationships below .05
${GCTA} \
--grm 1-GRM_FULL \
--thread-num 1 \
--grm-cutoff 0.05 \
--make-grm \
--out 1-GRM_FULL.RM5

## Done
echo `date` "2 - Calculate GRM - DONE" > ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 3 ## Calculate Principal Components ###################################
##########################################################################
if [ "$START_STEP" -le 3 ]; then
echo \### 3 - `date` \###
echo \### Calculate PCs \###
echo `date` "3 - Calculate PCs" >> ${UPDATE_FILE}

#####################################################
## Do PCA Analysis for Whole Genome
 # GCTA --pca
 ${GCTA} --pca 20 \
 --grm 1-GRM_FULL.RM5 \
 --thread-num 1 \
 --out 2-PCA_FULL

## Done
echo `date` "2 - Calculate GRM - DONE" > ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 4 ## Pull out Covariates ##############################################
##########################################################################
if [ "$START_STEP" -le 4 ]; then
echo \### 4 - `date` \###
echo \### Collect Covariates \###
echo `date` "4 - Collect Covariates" >> ${UPDATE_FILE}

## If Covariates are being Used
if [[ $USE_COVARS == TRUE ]]
then

#####################################################
## Compile Covariates into Correct Format
NEW_COV_PATH=${OUT_DIR}/${COVS_FILENAME}_FULL.txt
echo ${NEW_COV_PATH}
Rscript ${PULL_COVS} ${COVS_COMMAND} ${OUT_DIR}/2-PCA_FULL.eigenvec ${COV_PATH} ${NEW_COV_PATH}

fi

## Done
echo `date` "4 - Collect Covariates - DONE" > ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## 5 ## Estimate Heritability ############################################
##########################################################################
if [ "$START_STEP" -le 5 ]; then
echo \### 5 - `date` \###
echo \### Estimate Heritability \###
echo `date` "5 - Estimate Heritability" >> ${UPDATE_FILE}

## Loop through Phenotypes
for file in `cat ${PHENO_DIR}/${PHENO_FILE_LIST}`
do
EST_OUT=${OUT_DIR}/3-REML_file

## Run GCTA to get Heritability Estimates
 # If Covariates are Specified
if [[ $USE_COVARS == TRUE ]]
then
${GCTA} \
--grm 1-GRM_FULL.RM5 \
--pheno ${PHENO_DIR}/${file} \
--qcovar ${NEW_COV_PATH} \
--reml \
--reml-maxit 1000 \
--reml-est-fix \
--reml-pred-rand \
--out ${EST_OUT}
else
 # If Covariates are NOT Specified
${GCTA} \
--grm 1-GRM_FULL.RM5 \
--pheno ${PHENO_DIR}/${file} \
--reml \
--reml-maxit 1000 \
--reml-est-fix \
--reml-pred-rand \
--out ${EST_OUT}
fi

done # Close Phenotype Loop

## Done
echo `date` "5 - Estimate Heritability - DONE" > ${UPDATE_FILE}
printf "V\nV\nV\nV\nV\nV\nV\nV\n"
fi
##########################################################################
## END OF DOC ############################################################
##########################################################################












