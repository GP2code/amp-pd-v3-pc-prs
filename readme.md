# AMP-PD version 3: PC and PRS calculation

Please see **main.ipynb** for the actual processing of data. 

We do not share the individual level data here. Please visit https://amp-pd.org/ for the data access.

## PC calculation and population assignment
AMP-PD version 3 wgs genetics data of 10,418 samples were merged with the 1000 Genomes Project Phase 3 reference panel. Using the mean +/- 6SD threshold of PC1-5, the ancestry {EUR, EAS, AFR, OTHER} was assigned.

Ancestry assignment results:

* EUR 10119
* EAS 32
* AFR 111
* OTHER 156

![Population Plot](/ancestry/PopulationPlot_w_1kg.png)

## PRS calculation
The polygenic risk scores (PRS) were calculated based on the PD GWAS summary statistics from Nalls et al. 2019 (PMID: 31701892). The weight file was downloaded from  https://github.com/neurogenetics/genetic-risk-score. The file contains the weights of 90 risk SNPs but two SNPs were not available in the AMP-PD genetics data. We generated various versions of the PRSs and in some versions, the missing SNPs-rs11578699 was substituted by rs11577197 (perfect LD) and rs3742785 was substituted by rs10134885 (R2=0.989). The PRSs were normalized (N[0, 1]) across the AMP-PD for the ease of use. 

PRS variations:
* PRS88 -> Full (but missing 2 SNPs)
* PRS83 -> PRS88 Without GBA, LRRK2 loci
* PRSp90 -> Used proxy above
* PRSp85 -> Used proxy above without GBA, LRRK2 loci
* PRSp87 -> Used proxy above without GBA locus
* PRSp88 -> Used proxy above without LRRK2 locus


**Pleaes note that the PRS weights were derived from the European GWAS results. Therfore, the PRSs may be invalid for those with non-European ancestries.**