RUN_GCTA

## Names/Paths
echo \### Defining Set Variables and Paths at `date` \###
DATE=20150205_1M
HOME_DIR=/projects/janssen/Psych/GCTA
cd ${HOME_DIR}

## Files
VAR_FILE=1M_white_hispanic_scz_sca_qc
VAR_DIR=/projects/janssen/Psych/Data/Genotyped/
PHENO_DIR=/projects/janssen/Psych/Pheno
PHENO_FILE=Full_Table.txt
PHENO_NAME_LIST=PHENO_NAMES.txt
COV_FILE=Full_Table.txt
COVS=`echo BL_PANSS AGE_DIAG`
PC_COUNT=0
START_STEP=6

COVS=`echo "$COVS" | sed 's/ /QQQ/g'`

########################################
## Run The Script
/projects/janssen/Psych/Scripts/GCTA/GCTA.sh \
${DATE} \
${HOME_DIR} \
${VAR_FILE} \
${VAR_DIR} \
${PHENO_DIR} \
${PHENO_FILE} \
${PHENO_NAME_LIST} \
${COV_FILE} \
${COVS} \
${PC_COUNT} \
${START_STEP}

################################################################################
################################################################################
################################################################################

## Files
VAR_FILE=PsychChip_R092670-PSY-3006_R092670_arm

########################################
## Run The Script
/projects/janssen/Psych/Scripts/GCTA/GCTA.sh \
${DATE} \
${HOME_DIR} \
${VAR_FILE} \
${VAR_DIR} \
${PHENO_DIR} \
${PHENO_FILE} \
${PHENO_NAME_LIST} \
${COV_FILE} \
${COVS} \
${PC_COUNT} \
${START_STEP}