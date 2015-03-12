RUN_GCTA

#!/bin/bash
#PBS -N 20150310_GCTA
#PBS -q hotel
#PBS -v QOS=10
#PBS -l walltime=10:00:00
#PBS -m abe
#PBS -l nodes=1:ppn=1
#PBS -o 20150310_GCTA.run.oe
#PBS -j oe
#PBS -M k26stan@yahoo.com
#PBS -V
#PBS -A schork-group

echo "<startTime>"`date`"</startTime>"
echo "<output>"

## Names/Paths
echo \### Defining Set Variables and Paths at `date` \###
DATE=20150310_ALL5b
HOME_DIR=/projects/janssen/Heritability
cd ${HOME_DIR}

## Files
VAR_FILE=BED_MAF5.ALL
VAR_DIR=/projects/janssen/VCFs/PLINK
PHENO_DIR=/projects/janssen/ASSOCIATION/PH-PHENOTYPES
PHENO_FILE=20150310_Single_Pheno_Table.txt
PHENO_NAME_LIST_FILE=Manuscript_Phenos_Covs.txt
PHENO_NAME_LIST_DIR=/projects/janssen/Heritability
COV_FILE=20150310_Single_Pheno_Table.txt
PC_COUNT=0
START_STEP=1

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
${PHENO_NAME_LIST_FILE} \
${PHENO_NAME_LIST_DIR} \
${COV_FILE} \
${PC_COUNT} \
${START_STEP}





foo[@]:1:







#############################################################
#############################################################
#############################################################
#############################################################



















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