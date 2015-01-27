#!/bin/bash
GCTA=/gpfs/group/schork/nwineing/gcta/gcta64
DIR=/gpfs/group/schork/nwineing/jnj/gcta
PHENOARRAY[0]=CGIS-combined-european
PHENOARRAY[1]=CGIS-er_oros-european
PHENOARRAY[2]=CGIS-palmitate-european
PHENOARRAY[3]=nPANSS-combined-european
PHENOARRAY[4]=nPANSS-er_oros-european
PHENOARRAY[5]=nPANSS-palmitate-european
PHENOARRAY[6]=pPANSS-combined-european
PHENOARRAY[7]=pPANSS-er_oros-european
PHENOARRAY[8]=pPANSS-palmitate-european
PHENOARRAY[9]=PANSS-combined-european
PHENOARRAY[10]=PANSS-er_oros-european
PHENOARRAY[11]=PANSS-palmitate-european
for i in {0..11}
do
PHENO=${PHENOARRAY[${i}]}
#mkdir ${PHENO}
#${GCTA} --grm ${DIR}/combined --pca 10 --out ${DIR}/${PHENO}/${PHENO}
#R --no-save ${PHENO} < pcs.r
${GCTA} --grm ${DIR}/combined.info5 --pheno ${DIR}/datasets/${PHENO}.pheno --qcovar ${DIR}/${PHENO}/${PHENO}.covar --reml --reml-maxit 1000 --reml-est-fix --reml-pred-rand --out ${DIR}/${PHENO}/result.info5
done

