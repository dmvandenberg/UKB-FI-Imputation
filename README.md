# UKB-FI-Imputation
Here we list scripts used in the paper "Imputation of fluid intelligence scores reduces ascertainment bias and increases power for analyses of common and rare variants"
for a detailed description of analyses and decisions please refer to our paper.

## Imputation
Scripts regarding imputation

- **[prepare_FI_pheno.R](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/imputation/prepare_FI_pheno.R)**: is the script used to standardize FI measured and configure covariate files for GWAS with the appropriate values (age at time of measurement). This script requires UKB variables 34-0.0 and 52-0.0, and all instances of 20016, 20131, 21003, 20135. In addition it requires PC loadings for 25 ancestry PCs, and the sex and array used for genotyping for all individuals.
- **[FI_imputation.R](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/imputation/FI_imputation.R)**: is used for imputation. It leverages softImpute and requires access to various UKB datafield to impute FI. We have evaluated 2 sets of imputation variables for which the UKB variable ID's can be found [here](https://github.com/dmvandenberg/UKB-FI-Imputation/tree/main/imputation/imputation-variables).
- **gwas_by_subtraction.R**: #TOADD
- **[rescale_combine.R](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/imputation/rescale_combine.R)**: is used to mega-analyze measured and imputed FI values after standardizing them separately.
- **[combine.metal](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/imputation/combine.metal)**: is used to meta-analyze GWAS of measured and imputed FI.

## Genetic analysis

- **[gwas.job](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/genetic-analyses/gwas.job)**: is the script used to run a GWAS using [GCTA](https://yanglab.westlake.edu.cn/software/gcta/). This requires access to genotype data, a sparse GRM, and one of the .pheno files created by our imputation scripts along with accompanying .qcov file
- **[COGENT_meta.metal](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/genetic-analyses/COGENT_meta.metal)**: is the [metal](https://genome.sph.umich.edu/wiki/METAL_Documentation) meta-analysis script used to meta-analyze the GWAS of COGENT and Combined FIS
  
### Common variant analyes

Genetic correlations and SNP heritabilities were computed using [LDSC](https://github.com/bulik/ldsc). The scripts we used for this are:

- **[munge.job](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/genetic-analyses/munge.job)**: Script to format following LDSC standards
- **[h2_rg.job](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/genetic-analyses/h2_rg.job)**: To compute SNP heritability and genetic correlations with reference GWAS

### Rare variant analyses

1. Quality control of rare genetic variants
   We used applets designed by Eugene Gardner to process the quality control of rare variants. For details, please refer to https://github.com/mrcepid-rap. specifically, we used mrc-splitbcf, mrc-filterbcf, mrc-makebgen and mrc-collapsevariants to obtain the variant   files for regenie. We filtered the variants by minor allele frequency (MAF), retaining only those with a MAF lower than 0.001%. We also annotated two additional variant effect predictor MPC and alphamissense to the missense variants using in-house written script [annotate-score.sh](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/rare-variant-analyses/annotate-score.sh).

   For PTVs, only high-confidence PTVs defined by LOFTEE were retained. For missense variants, a damaging missense variant set was created by including variants (Alphamissense >= 0.56, & REVEL >= 0.5 & MPC >=2).  

2. format the variant file for regenie
   We specifically use <output_prefix>..STAAR.variants_table.tsv from the Output Tarball of mrc-collapsevariants. We used [UKB_fileformat.R](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/rare-variant-analyses/UKB_fileformat.R) to produce the setlist and annotation files required to run gene-based testing in Regenie.

3. Running regenie on UKB research analysis platform (UKB-RAP)
   After files were prepared for regenie, we used two steps to run gene-based testing in Regenie on UKB-RAP. We used swiss-army-knife to run regenie.
   a) Pre-filtering before running regenie.

   Merge common genetic variant files
   `dx run swiss-army-knife -iin="default.txt" -icmd='cp /project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22418_c[1-9]* . ;ls *.bgen | sed -e 's/.bgen//g'> files_to_merge.txt;plink --merge-list files_to_merge.txt --make-bed --autosome-xy --out ukb22828_c1_22_v2_merged;rm files_to_merge.txt;' --tag='Step1' --instance-type "mem1_ssd1_v2_x16" --destination='/project/Data/' --brief --yes`

   Quality control the SNPs
   `dx run swiss-army-knife -iin="/project/Data/ukb22418_c1_22_v2_merged.bed" -iin="/project/Data/ukb22418_c1_22_v2_merged.bim" -iin="/project/Data/ukb22418_c1_22_v2_merged.fam" -icmd="plink2 --bfile ukb22418_c1_22_v2_merged --autosome --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1 --write-snplist --write-samples --no-id-header --out qc_pass" --tag="Array QC" --instance-type "mem2_ssd2_v2_x96" --destination="/project/Data/" --brief --yes`

   b) Example: regenie step1

   Input files:
   - genetic data file: bed/bim/fam file
   - phenotype data file: ukb_phenotype.txt including focal phenotype and covariates
   - covariate data file: same with phenotype data file
   - optional: list of variants to keep (qc_pass.snplist)
   - output file: Predictions files (.loco.gz)
  
   `dx run swiss-army-knife -iin="ukb22418_c1_22_v2_merged.bed" -iin="ukb22418_c1_22_v2_merged.bim" -iin="ukb22418_c1_22_v2_merged.fam" -iin="qc_pass.snplist" -iin="ukb_phenotype.txt" -icmd="regenie --step 1 --bed ukb22418_c1_22_v2_merged --phenoFile ukb_phenotype.txt --covarFile ukb_phenotype.txt --extract qc_pass.snplist --phenoColList fluid.intel --covarColList age,birth.year,pca{1:25},age.squared,agebysex,age.squaredbysex,sex --catCovarList center --maxCatLevels 30 --lowmem --bsize 1000 --qt --apply-rint --gz --verbose --threads 16 --out fis_step1" --instance-type "mem1_ssd1_v2_x16" --destination='/result/' --brief --yes`

    c) Example: regenie step2 (eg.chr1)

   Input files:
   - genetic data file: BGEN file and Sample file corresponding to input BGEN file
   - phenotype data file: ukb_phenotype.txt including focal phenotype and covariates
   - covariate data file: same with phenotype data file
   - prediction files: ...pred.list and ...loco.gz from step1 results
   - annotation input files: ...setlist.txt and ...annotations.txt (see above for file preparations)
   - mask files: mask.txt specifies which annotation categories should be combined into masks
   - sample files: ..BOLT.sample (BOLT-ready bgen sample file of per-gene 'genotypes' from mrc-collapsevariants
   - output file: regenie results files
  
   `dx run swiss-army-knife -iin="/filtered-bgen/chr1.filtered.bgen" -iin="/project/Data/HC-PTV-chr1.22.BOLT.sample"   -iin="/project/Data/ukb_phenotype.txt" -iin="/project/Data/fis_step1_pred.list" -iin="/project/Data/fis_step1_1.loco.gz" -iin="/project/step2_lof/chr1_setlist.txt" -iin="/project/step2_lof/chr1_annotations.txt" -iin="/project/step2_lof/mask.txt" -icmd="regenie --step 2 --bgen chr1.filtered.bgen --sample HC-PTV-chr1.22.BOLT.sample --phenoFile ukb_phenotype.txt --covarFile ukb_phenotype.txt --phenoColList fluid.intel --covarColList age,birth.year,pca{1:25},age.squared,agebysex,age.squaredbysex,sex --catCovarList center --maxCatLevels 30 --firth --approx --pred fis_step1_pred.list  --anno-file chr1_annotations.txt --set-list chr1_setlist.txt --mask-def mask.txt --aaf-bins 1 --vc-tests skato-acat,acato-full --bsize 200 --qt --apply-rint --threads 2 --gz --verbose --out chr1_ptv" --tag='Step1' --instance-type "mem1_ssd1_v2_x16" --destination='/project/Data/' --brief --yes`
