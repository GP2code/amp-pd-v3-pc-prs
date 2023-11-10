#!bin/bash

# Take the plink1.9 binary file and calculate PCA

# Input variables
INPUT=$1
NTHREAD=$2

# If NTHREAD is empty, set it 2
if [ -z "$NTHREAD" ]
then
  NTHREAD=2
fi

## get prune.in
plink\
  --bfile ${INPUT}\
  --indep-pairwise 50 5 0.01\
  --out ${INPUT}_pruning\
  --threads ${NTHREAD}

## prune
plink\
  --bfile ${INPUT}\
  --extract ${INPUT}_pruning.prune.in\
  --make-bed\
  --out ${INPUT}_pruned\
  --threads ${NTHREAD}

## pca
plink\
  --bfile ${INPUT}_pruned\
  --out ${INPUT}_pca\
  --pca\
  --threads ${NTHREAD}
