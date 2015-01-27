#!/bin/bash
GCTA=/gpfs/group/schork/nwineing/gcta/gcta64
DIR=/gpfs/group/schork/nwineing/jnj/gcta
${GCTA} --bfile /gpfs/group/schork/nwineing/jnj/gcta/datasets/combined.info5 --autosome --make-grm --out ${DIR}/combined.info5

