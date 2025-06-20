#!/sw/arch/RHEL8/EB_production/2022/software/R/4.2.1-foss-2022a/bin/Rscript
#SBATCH --output=/home/dvandenberg/gwas-by-sub.log
#SBATCH --time=04:00:00
#SBATCH --partition=rome
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --nodes=1

# GenmicSEM Cholesky model for IQ and ImputedIQ GWAS
.libPaths("./GenomicSEM/rpackages")

library(GenomicSEM)
library(data.table)

# munge sumstats.
munge("./SavageJansen_2018_intelligence_metaanalysis.txt", 
      "./w_hm3.snplist",
      trait.names="CTG_FIS", 
      N=269867, 
      info.filter = 0.9, 
      maf.filter = 0.01)

munge("./imputed_FIS.fastGWA", 
      "./w_hm3.snplist",
      trait.names="IMP_FIS", 
      N=309128,
      info.filter = 0.9,
      maf.filter = 0.01
      )

# Produce S&V matrix
traits <- c("IMP_FIS.sumstats.gz","CTG_FIS.sumstats.gz")
sample.prev <- c(NA,NA)
population.prev <- c(NA,NA)
ld<-"./software/ldsc/eur_w_ld_chr"
wld <- "./software/ldsc/eur_w_ld_chr"
trait.names<-c("IMP_FIS", "CTG_FIS")

LDSCoutput <- ldsc(traits, 
                   sample.prev, 
                   population.prev, 
                   ld, 
                   wld, 
                   trait.names)

save(LDSCoutput, file="./LDSCoutputFISimpFIS.RData")

#-------------------------------------------------------------------
# Read in the summary statistics file
sumstats <- fread("SavageJansen_2018_intelligence_metaanalysis_hm3.txt", sep="\t")
print("file read")

# Ensure the columns are numeric
sumstats$BETA <- as.numeric(sumstats$BETA)
sumstats$SE <- as.numeric(sumstats$SE)
sumstats$P <- as.numeric(sumstats$P)
sumstats$N <- as.numeric(sumstats$N)
print("conversion succeeded")

# Remove rows with NA values in these columns
print(dim(sumstats[is.na(sumstats$BETA) | is.na(sumstats$SE) | is.na(sumstats$P) | is.na(sumstats$N),]))
sumstats <- sumstats[!is.na(sumstats$BETA) & !is.na(sumstats$SE) & !is.na(sumstats$P) & !is.na(sumstats$N),]

# Write the cleaned data back to a file
print("writing file")
fwrite(sumstats, "SavageJansen_2018_intelligence_metaanalysis_hm3.cleaned.txt", sep="\t", 
       quote=F, col.names=T, row.names=F)

# Clean SNPs for GWAS
files = c("SavageJansen_2018_intelligence_metaanalysis_hm3.cleaned.txt",'imputed_FIS1.v1.fastGWA')
ref = "reference.1000G.maf.0.005.txt"
trait.names<-c("CTG_FIS", "IMP_FIS")
se.logit = c(F,F)
info.filter = 0.6
maf.filter = 0.01

p_sumstats<-sumstats_dirk(files, ref, trait.names, se.logit, info.filter, maf.filter, OLS=c(F,F), linprob=NULL, N=c(269867,309102), betas=NULL)

save(p_sumstats, file="Sumstats.RData")


