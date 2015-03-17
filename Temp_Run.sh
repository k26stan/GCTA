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
DATE=20150316_Single_MAF5_ALL
HOME_DIR=/projects/janssen/Heritability
cd ${HOME_DIR}

## Files
VAR_FILE=BED_MAF5.ALL # FULL_BED_ALL/BED_MAF5.ALL
VAR_DIR=/projects/janssen/VCFs/PLINK
PHENO_DIR=/projects/janssen/ASSOCIATION/PH-PHENOTYPES
PHENO_FILE=20150310_Single_Pheno_Table.txt
PHENO_NAME_LIST_FILE=Manu_PhenoCovs_Single.txt
PHENO_NAME_LIST_DIR=/projects/janssen/Heritability
COV_FILE=20150310_Single_Pheno_Table.txt
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
