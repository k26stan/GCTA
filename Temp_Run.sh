## Notes for Updated GCTA Pipeline
 # Using new version of GCTA
# Using GREML-LDMS Approach
  # That is Robust to LD and MAF


## Names/Paths
echo \### Defining Set Variables and Paths at `date` \###
DATE=20150924_Test_Run
HOME_DIR=/projects/janssen/Heritability
cd ${HOME_DIR}

## Files
VAR_FILE=BED_FULL.ALL # PsychChip_R092670-PSY-3006_R092670_arm
VAR_DIR=/projects/janssen/VCFs/PLINK
PHENO_DIR=/projects/janssen/ASSOCIATION/PH-PHENOTYPES
PHENO_FILE=20150506_Derived_Pheno_Table.txt
PHENO_NAME_LIST_FILE=Manu_PhenoCovs_Derived.txt
PHENO_NAME_LIST_DIR=/projects/janssen/Heritability
COV_FILE=20150506_Derived_Pheno_Table.txt
PC_COUNT=0
START_STEP=3

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



############################################################################
############################################################################
## PSYCH COHORTS ###########################################################
############################################################################
############################################################################

########################################################
## 1M Cohort ##

## Names/Paths
echo \### Defining Set Variables and Paths at `date` \###
DATE=20150319_1M
HOME_DIR=/projects/janssen/Psych/GCTA
cd ${HOME_DIR}

## Files
VAR_FILE=1M_white_hispanic_scz_sca_qc # PsychChip_R092670-PSY-3006_R092670_arm
VAR_DIR=/projects/janssen/Psych/Data/Genotyped
PHENO_DIR=/projects/janssen/Psych/Pheno
PHENO_FILE=Full_Table.txt
PHENO_NAME_LIST_FILE=PhenoCovs_Table_NoSex.txt
PHENO_NAME_LIST_DIR=/projects/janssen/Psych/Pheno
COV_FILE=Full_Table.txt
PC_COUNT=0
START_STEP=4

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


########################################################
## PsychChip COHORT ##
DATE=20150319_Psych
VAR_FILE=PsychChip_R092670-PSY-3006_R092670_arm
PHENO_NAME_LIST_FILE=PhenoCovs_Table_NoSex.txt

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

############################################################################
############################################################################
## RA COHORT ###############################################################
############################################################################
############################################################################

########################################################
## DERIVED PHENOTYPES ##
MAF5/MAF1
ALL/SNP

## Names/Paths
echo \### Defining Set Variables and Paths at `date` \###
DATE=20150316_Derived_MAF5_ALL
HOME_DIR=/projects/janssen/Heritability
cd ${HOME_DIR}

## Files
VAR_FILE=BED_MAF5.ALL # FULL_BED_ALL/BED_MAF5.ALL
VAR_DIR=/projects/janssen/VCFs/PLINK
PHENO_DIR=/projects/janssen/ASSOCIATION/PH-PHENOTYPES
PHENO_FILE=20150313_Derived_Pheno_Table.txt
PHENO_NAME_LIST_FILE=Manu_PhenoCovs_Derived.txt
PHENO_NAME_LIST_DIR=/projects/janssen/Heritability
COV_FILE=20150313_Derived_Pheno_Table.txt
PC_COUNT=0
START_STEP=1

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

########################################################
## SINGLE PHENOTYPES ##
MAF5/MAF1
ALL/SNP

## Names/Paths
echo \### Defining Set Variables and Paths at `date` \###
DATE=20150316_Derived_MAF1_ALL
HOME_DIR=/projects/janssen/Heritability
cd ${HOME_DIR}

## Files
VAR_FILE=BED_MAF1.ALL # FULL_BED_SNP/BED_MAF5.SNP # 
VAR_DIR=/projects/janssen/VCFs/PLINK
PHENO_DIR=/projects/janssen/ASSOCIATION/PH-PHENOTYPES
PHENO_FILE=20150313_Derived_Pheno_Table.txt
PHENO_NAME_LIST_FILE=Manu_PhenoCovs_Derived.txt
PHENO_NAME_LIST_DIR=/projects/janssen/Heritability
COV_FILE=20150313_Derived_Pheno_Table.txt
PC_COUNT=0
START_STEP=7

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
