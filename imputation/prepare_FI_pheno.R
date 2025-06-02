#!/sw/arch/RHEL8/EB_production/2022/software/R/4.2.1-foss-2022a/bin/Rscript
#SBATCH --output=./prepare_fis.log
#SBATCH --time=01:00:00
#SBATCH --partition=rome
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --nodes=1

# load packages 
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(lubridate))

# ---- 0. SETTING SOME PARAMETERS ---- #
# Extracted UKB variables (all instances) 20016, 20131, 21003, 20135, 34-0.0 and 52-0.0
raw_FIS_with_ages <- "/home/dvandenberg/project_folder/iq-from-scratch/fis_measures/selected_phenos.csv"

# Columns: EID and loadings for 25 genectic PCs
EU_genetic_PCs_file <- "./age.age_sq.25PCs.qcov"

# Columns: Number of y chromosomes and array used for genotyping
sex_array_file <- "./sex.array.cov"

# Output folder
output_folder <- "./FI_phenos_prepared/"

# SELECT A FIS APPROACH:
#  - "earliest": select earliest FIS measure
#  - "average": select average FIS measure
#  - "all5": select all 5 FIS measures
#  - "FIS1": select ONLY FIS1 measure

FIS_approach <- "average"
residualize_bool <- TRUE # Set to FALSE if you don't want to residualize

# ---- 1. READ AND FILTER DATA ---- #

# Read FIS phenotypes
raw_FIS_with_ages <- fread(raw_FIS_with_ages, header=TRUE, sep=",", quote='"')

# Parse genetic PC file
EU_ancestry <- read.table(EU_genetic_PCs_file, header=FALSE)
colnames(EU_ancestry) <- c("eid", "fid", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
                           "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", 
                           "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", 
                           "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", 
                           "PC25")
EU_ancestry$'fid' <- NULL

# Parse sex, array file
sex.array <- read.table(sex_array_file, header=FALSE)
colnames(sex.array) <- c("eid", "fid", "sex", "array")
sex.array$'fid' <- NULL

# Combine everything into one dataframe
covariates <- merge(EU_ancestry, sex.array, by="eid", all = FALSE)
FIS_with_ages <- merge(raw_FIS_with_ages, covariates, by="eid", all = FALSE)

# Remove separate dataframes to save space
rm(EU_ancestry)
rm(sex.array)
rm(covariates)
rm(EU_genetic_PCs_file)
rm(wei_download)

#### Here we have 455665 genetically EU individuals

# ---- 2. PROCESS AGE DATA ---- #
# Assesment visits have an age at visit datafield, online questionaire only has date taken (20135)
# 2.1 Convert UNIX timestamp to date for datafields 20135
FIS_with_ages$`20135-0.0` <- as.POSIXct(FIS_with_ages$`20135-0.0`, origin="1970-01-01", tz="UTC")
FIS_with_ages$`20135-1.0` <- as.POSIXct(FIS_with_ages$`20135-1.0`, origin="1970-01-01", tz="UTC")

# 2.2 Combine birth year and month into birth month column
FIS_with_ages$`34-0.0` <- as.integer(FIS_with_ages$`34-0.0`)
FIS_with_ages$`52-0.0` <- as.integer(FIS_with_ages$`52-0.0`)
FIS_with_ages$birth_month <- paste(sprintf("%04d", FIS_with_ages$`34-0.0`), 
                                   sprintf("%02d", FIS_with_ages$`52-0.0`), sep = "-")
FIS_with_ages$`34-0.0` <- NULL
FIS_with_ages$`52-0.0` <- NULL

# 2.3 Compute approximate age at measurement 20135
#  by calculating difference between birth month and measurement date
FIS_with_ages$`20135-0.0` <- as.Date(FIS_with_ages$`20135-0.0`, format="%Y-%m-%d")
FIS_with_ages$`20135-1.0` <- as.Date(FIS_with_ages$`20135-1.0`, format="%Y-%m-%d")
FIS_with_ages$`birth_month` <- as.Date(paste0(FIS_with_ages$`birth_month`,"-01"), format="%Y-%m-%d")
FIS_with_ages$`20135-0.0` <- as.numeric(round(interval(FIS_with_ages$birth_month, FIS_with_ages$`20135-0.0`) / years(1)))
FIS_with_ages$`20135-1.0` <- as.numeric(round(interval(FIS_with_ages$birth_month, FIS_with_ages$`20135-1.0`) / years(1)))
FIS_with_ages$`birth_month` <- NULL

# ---- 3. PROCESS FIS DATA ---- #
# 1) truncate FI scores of 14
# 2) Residualize FIS if required
#   We perform the following regression:
#   FIS + intercept ~ age + age^2

# Make sex factor as they are categorical
FIS_with_ages$sex <- factor(FIS_with_ages$sex)

# First truncate a FIS score of 14 to 13
FIS_with_ages$`20191-0.0`[FIS_with_ages$`20191-0.0` > 13] <- 13
FIS_with_ages$`20191-1.0`[FIS_with_ages$`20191-1.0` > 13] <- 13

# Declare residualization function
residualize <- function(data, fis, age) {
  # make sure labels are being interpreted correctly
  fis <- paste0("`",fis,"`")
  age <- paste0("`",age,"`")
  # Create formula
  formula <- as.formula(paste0(fis, " ~ ", age, " + I(", age, "^2)"))
  print(formula)
  # Run model
  model <- lm(formula, data = data)
  residuals <- resid(model)
  # Add model intercept back into residuals to allow for different means
  intercept <- coef(model)[1]
  print(paste("Intercept:", intercept))
  residuals <- residuals + intercept
  # return residuals and indices for which row they are
  return(list(residuals = residuals, indices = as.integer(names(residuals)), model=model))
}

if (residualize_bool == TRUE){
  # Residualize proportions
  resid_FIS1 <- residualize(FIS_with_ages, "20016-0.0", "21003-0.0")
  resid_FIS2 <- residualize(FIS_with_ages, "20016-1.0", "21003-1.0")
  resid_FIS3 <- residualize(FIS_with_ages, "20016-2.0", "21003-2.0")
  resid_FIS4 <- residualize(FIS_with_ages, "20191-0.0", "20135-0.0")
  resid_FIS5 <- residualize(FIS_with_ages, "20191-1.0", "20135-1.0")
  
  # Get residuals back into data frame
  FIS_with_ages$`20016-0.0` <- NA_real_
  FIS_with_ages$`20016-0.0`[resid_FIS1$indices] <- resid_FIS1$residuals 
  FIS_with_ages$`20016-1.0` <- NA_real_
  FIS_with_ages$`20016-1.0`[resid_FIS2$indices] <- resid_FIS2$residuals
  FIS_with_ages$`20016-2.0` <- NA_real_
  FIS_with_ages$`20016-2.0`[resid_FIS3$indices] <- resid_FIS3$residuals
  FIS_with_ages$`20191-0.0` <- NA_real_
  FIS_with_ages$`20191-0.0`[resid_FIS4$indices] <- resid_FIS4$residuals
  FIS_with_ages$`20191-1.0` <- NA_real_
  FIS_with_ages$`20191-1.0`[resid_FIS5$indices] <- resid_FIS5$residuals
}

# ---- 4. CLEAN AND SPLIT UP DATAFRAME ---- #
# We want to end up with:
#  1) a set of (residualized) FIS phenotypes to use
#  2) An accompanying .qcov file
#  3) A sex.array.cov file (the input file is good as is)

# set up some subsets
FIS <- FIS_with_ages[, c("eid", "20016-0.0", "20016-1.0", "20016-2.0", "20191-0.0", "20191-1.0")]
ages <- FIS_with_ages[, c("eid", "21003-0.0", "21003-1.0", "21003-2.0", "20135-0.0", "20135-1.0")]
base_covariate <- FIS_with_ages[, c("eid", "eid", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6",
                                    "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", 
                                    "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", 
                                    "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", 
                                    "PC25")]
colnames(base_covariate)[2] <- "fid"

# Create FIS phenos and accompanying covariate file
if (FIS_approach == "FIS1" | FIS_approach == "all5"){
  # make qcov (FIS1 and all5 have the same covariate files)
  qcov <- merge(base_covariate, ages[, c("eid", "21003-0.0")], by="eid", all.x = TRUE)
  # Compute age^2, age*sex, age^2*sex
  qcov$'age.sq' <- I(qcov$"21003-0.0"^2)
  qcov$'age.sex' <- I(qcov$'21003-0.0'*(as.numeric(FIS_with_ages$sex)-1))
  qcov$'age.sq.sex' <- I(qcov$'age.sq'*(as.numeric(FIS_with_ages$sex)-1))
  write.table(qcov, paste0(output_folder, "FIS1.qcov"), col.names = FALSE, row.names = FALSE)
  
  if (FIS_approach == "FIS1") {
    # make FIS
    FIS <- FIS[, c("eid", "20016-0.0")]
    write.table(FIS, paste0(output_folder, "FIS1.pheno"), sep="\t", row.names = FALSE)
  } else {
    write.table(FIS, paste0(output_folder, "all5FIS.pheno"), sep="\t", row.names = FALSE)
  }
}

if (FIS_approach == "average"){
  # Average ages
  avg_ages <- ages
  avg_ages$mean_age <- avg_ages %>% select(-eid) %>% rowMeans(na.rm = TRUE)
  # Add average age to base covariates
  avg_qcov <- merge(base_covariate, avg_ages[, c("eid", "mean_age")], by="eid", all.x = TRUE)
  # Compute age^2, age*sex, age^2*sex
  avg_qcov$'age.sq' <- I(avg_qcov$"mean_age"^2)
  avg_qcov$'age.sex' <- I(avg_qcov$"mean_age"*(as.numeric(FIS_with_ages$sex)-1))
  avg_qcov$'age.sq.sex' <- I(avg_qcov$'age.sq'*(as.numeric(FIS_with_ages$sex)-1))
  write.table(avg_qcov, paste0(output_folder, "average_FIS.qcov"), col.names = FALSE, row.names = FALSE)
  
  # Compute average FIS
  avg_FIS <- FIS
  avg_FIS$mean_FIS <- avg_FIS %>% select(-eid) %>% rowMeans(na.rm = TRUE)
  write.table(avg_FIS[, c("eid", "mean_FIS")], paste0(output_folder, "average_FIS.pheno"), sep="\t", row.names = FALSE)
}

if (FIS_approach == "earliest"){
  earliest_FIS <- FIS
  earliest_measures <- t(apply(earliest_FIS[, -which(names(earliest_FIS) == "eid")], 1, function(x) {
    idx <- which(!is.na(x))[1]
    if (!is.na(idx)) {
      c(value = x[idx], index = idx+1)
    } else {
      c(value = NA, index = 2)
    }
  }))
  colnames(earliest_measures) <- c("value", "index")
  earliest_FIS$earliest <- earliest_measures[, "value"]
  
  earliest_age <- ages
  earliest_age$earliest <- earliest_measures[, "index"]
  earliest_age$earliest <- apply(earliest_age, 1, function(x) {
    age_column_index = as.numeric(x["earliest"])
    x[age_column_index]
  })
  
  write.table(earliest_FIS[, c("eid", "earliest")], paste0(output_folder, "earliest_FIS.pheno"), sep="\t", row.names = FALSE)
  
  earliest_qcov <- merge(base_covariate, earliest_age[, c("eid", "earliest")], by="eid", all.x = TRUE)
  earliest_qcov$'age.sq' <- I(earliest_qcov$"earliest"^2)
  earliest_qcov$'age.sex' <- I(earliest_qcov$"earliest"*(as.numeric(FIS_with_ages$sex)-1))
  earliest_qcov$'age.sq.sex' <- I(earliest_qcov$'age.sq'*(as.numeric(FIS_with_ages$sex)-1))
  write.table(earliest_qcov, paste0(output_folder, "earliest_FIS.qcov"), sep="\t", row.names = FALSE, col.names = FALSE)
}