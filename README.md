# UKB-FI-Imputation
Here we list scripts used in the paper "Imputation of fluid intelligence scores reduces ascertainment bias and increases power for analyses of common and rare variants"
for a detailed description of analyses and decisions please refer to our paper.

## A global overview of scripts in the general order in which to use them:

### Imputation
All scripts regarding imputation

- **prepare_FI_pheno.R**: is the script used to standardize FI measured and configure covariate files for GWAS with the appropriate values (age at time of measurement). This script requires UKB variables 34-0.0 and 52-0.0, and all instances of 20016, 20131, 21003, 20135. In addition it requires PC loadings for 25 ancestry PCs, and the sex and array used for genotyping for all individuals.
- **FI_imputation.R**: is used for imputation. It leverages softImpute and requires access to various UKB datafield to impute FI. We have evaluated 2 sets of imputation variables for which the UKB variable ID's can be found in the folder `imputation/imputation-variables`.
- **gwas_by_subtraction.R**: #TOADD
- **rescale_combine.R**: is used to mega-analyze measured and imputed FI values after standardizing them separately.
- **combine.metal**: is used to meta-analyze GWAS of measured and imputed FI.

### Genetic analysis

- **gwas.job**: is the script used to run a GWAS using [GCTA](https://yanglab.westlake.edu.cn/software/gcta/). This requires access to genotype data, a sparse GRM, and one of the .pheno files created by our imputation scripts along with accompanying .qcov file
- **COGENT_meta.metal**: is the [metal](https://genome.sph.umich.edu/wiki/METAL_Documentation) meta-analysis script used to meta-analyze the GWAS of COGENT and Combined FIS
  
#### Common variant analyes

Genetic correlations and SNP heritabilities were computed using [LDSC](https://github.com/bulik/ldsc). The scripts we used for this are:

- **munge.job**: Script to format following LDSC standards
- **h2_rg.job**: To compute SNP heritability and genetic correlations with reference GWAS

#### Rare variant analyses

- #TOADD


