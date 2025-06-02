#!/sw/arch/RHEL8/EB_production/2022/software/R/4.2.1-foss-2022a/bin/Rscript
#SBATCH --output=
#SBATCH --time=01:00:00
#SBATCH --partition=rome
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --nodes=1

suppressMessages(library(dplyr))

# ------- FILES ------- #
MEASURED <- "./FI-imputation/measured_average_FIS.pheno"
IMPUTED <- "./FI-imputation/imputed_average_FIS.pheno"
AVG_QCOV <- "./FI_phenos_prepared/average_FIS.qcov"
IMP_QCOV <- "./FI_phenos_prepared/FIS1.qcov"

# ------ COMBINE ------ #
# Read data
MEASURED <- read.table(MEASURED, sep="\t", header=TRUE)
message(paste0(length(MEASURED$EID), " individuals in measured file"))
IMPUTED <- read.table(IMPUTED, sep="\t", header=TRUE)
message(paste0(length(IMPUTED$EID), " individuals in imputed file"))
message(paste0(length(intersect(IMPUTED$EID, MEASURED$EID)), " individuals in both files (should be 0)"))

# Scale separately
MEASURED$PHENO <- scale(MEASURED$PHENO)
IMPUTED$PHENO <- scale(IMPUTED$PHENO)

# Combine
MEGA <- bind_rows(MEASURED, IMPUTED)
message(paste0(length(MEGA$EID), " individuals in MEGA file"))

# ---- COMBINE QCOV ---- #
AVG_QCOV <- read.table(AVG_QCOV, header=FALSE)
IMP_QCOV <- read.table(IMP_QCOV, header=FALSE)

names(AVG_QCOV)

UNIQUE_IMP_QCOV <- IMP_QCOV %>%
  filter(!V1 %in% AVG_QCOV$V1)

COMBINED_QCOV <- bind_rows(AVG_QCOV, UNIQUE_IMP_QCOV)
head(COMBINED_QCOV, 10)
message(paste0(length(COMBINED_QCOV$V1), " individuals in combined QCOV file"))

# ------- WRITE ------- #
write.table(MEGA, "average_mega.pheno", row.names = FALSE, sep='\t')
write.table(COMBINED_QCOV, "average_mega.qcov", row.names = FALSE, sep='\t', col.names = FALSE)