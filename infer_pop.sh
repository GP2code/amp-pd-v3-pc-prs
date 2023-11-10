#!bin/bash

# Take the plink2 (hg38) file
# Reduce to good quality biallelic SNPs in autosomes of maf 0.05, geno 0.05, mac 20
# Reduce to good quality sample with mind 0.05
# Standardize to hg38.fa.gz and rename
# Population reference is from 1000KG phase 3. 
## the reference was processed with https://github.com/hirotaka-i/1kg_ref


# Input variables
INPUT=$1
REF_PATH=$2
NTHREAD=$3
WORKDIR=$4

# If NTHREAD is empty, set it 2
if [ -z "$NTHREAD" ]
then
  NTHREAD=2
fi

# if WORKDIR is empty, set to current directory
if [ -z "$WORKDIR" ]
then
  WORKDIR=$(pwd)
fi

# VARS
FILE=$(basename ${INPUT})
REF_FILE=$(basename ${REF_PATH})
CURRENT_DIR=$(pwd)

# hg38 fasta file
wget ${WORKDIR} https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# reduce to good quality biallelic SNPs in autosomes
plink2\
  --autosome\
  --fa ${WORKDIR}/hg38.fa.gz\
  --geno 0.05\
  --mac 20\
  --maf 0.003\
  --make-pgen\
  --max-alleles 2\
  --mind 0.05\
  --normalize\
  --out ${WORKDIR}/${INPUT}_temp\
  --pfile ${INPUT}\
  --ref-from-fa force\
  --snps-only just-acgt\
  --sort-vars\
  --threads ${NTHREAD}

# move to the WORKDIR
cd ${WORKDIR}

# rename to reflect current chr:pos:ref:alt
plink2\
  --make-bed\
  --out ${INPUT}_snp\
  --pfile ${INPUT}_temp\
  --set-all-var-ids 'chr@:#:$r:$a'\
  --threads ${NTHREAD}

# merge with reference
## cp 1kg ref
cp ${REF_PATH}.fam ${REF_FILE}.fam
cp ${REF_PATH}.bim ${REF_FILE}.bim
cp ${REF_PATH}.bed ${REF_FILE}.bed

# premerge
plink\
 --bfile ${INPUT}_snp\
 --extract all_hg38_filtered_chrpos.bim\
 --make-bed\
 --out ${INPUT}_snp_premerge\
  --threads ${NTHREAD}

## need to use plink1.9 for this operation because plink2 does not support
## non-concatenating merging at this point
plink\
  --bfile ${INPUT}_snp_premerge\
  --bmerge ${REF_FILE}\
  --geno 0.01\
  --maf 0.01\
  --make-bed\
  --out ${INPUT}_snp_ref\
  --threads ${NTHREAD}

# prune_pca (plink1.9)
bash ../prune_pca.sh ${INPUT}_snp_ref ${NTHREAD}

# Save the analysis results and figures
python3 ../process_pca_results.py --input ${INPUT}_snp_ref_pca

# Move the file to the results directory
mkdir -p $CURRENT_DIR/ancestry
cp genetic_ancestry_* *.png $CURRENT_DIR/ancrestry