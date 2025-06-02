#!/sw/arch/RHEL8/EB_production/2022/software/R/4.2.1-foss-2022a/bin/Rscript
#SBATCH --output=./FI-imputation/imputation_job.log
#SBATCH --time=01:00:00
#SBATCH --partition=rome
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --nodes=1

# load packages 
suppressMessages(library(data.table))
suppressMessages(library(softImpute))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))

# ---- 0. SETTING PARAMETERS ---- #
ukb_phenotype_file <- "./ukb.phenotypes.csv" 
selected_phenotypes_file <- "./final_phenotype_selection.csv"
prepared_FIS_pheno_file <- "./FI-phenos-prepared/average_FIS.pheno"
output_folder <- "./FI-imputation/"

message("Selected parameters:")
message(paste0("prepared FIS file -- ", prepared_FIS_pheno_file))
message(paste0("phenotype selection -- ", selected_phenotypes_file))
message(paste0("output destination -- ", output_folder))

# Set random seed for reproducability
set.seed(1001)

# ---- 1. READING IN THE DATA ---- #
# read prepared FIS phenos (scaled, corrected and subsetted)
FIS_measures <- read.table(prepared_FIS_pheno_file, header=TRUE, sep="\t", check.names=FALSE)
included_FIS <- colnames(FIS_measures)[-1]

message("1. Loading data")
# Read file containing the ids of the phenotypes for imputation and then extract these from the UKB phenotype file
selected_phenotypes <- read.csv(selected_phenotypes_file, header=T)
selected_phenotypes <- as.vector(selected_phenotypes$pheno[selected_phenotypes$pheno!=""])

# Read UKB phenotype file headers and select only imputation phenotypes by id
UKBB.raw.phenotype.headers <- fread(ukb_phenotype_file, quote='"', sep=',', nrows = 0)
pheno_indices <- which(names(UKBB.raw.phenotype.headers) %in% selected_phenotypes)
message(paste0(" ... found ",length(selected_phenotypes), " phenotypes"))

UKBB.phenotypes <- fread(ukb_phenotype_file, quote='"', sep=',',select=pheno_indices)
dim(UKBB.phenotypes)

# Remove obsolete data
rm(selected_phenotypes)
rm(UKBB.raw.phenotype.headers)
rm(pheno_indices)

save.image(file = paste0(output_folder, "imputation_workspace_1.RData"))

# ---- 2. PROCESSING THE IMPUTATION FIELDS ---- #
message("2. Cleaning phenotypes")

# Perform phenotype cleaning
UKBB.phenotypes.temp <- UKBB.phenotypes
UKBB.phenotypes.temp[UKBB.phenotypes.temp==-7] <- 100
UKBB.phenotypes.temp1 <- UKBB.phenotypes.temp[,-c("5084-0.0","5085-0.0","5201-0.0")]
UKBB.phenotypes.temp2 <- UKBB.phenotypes.temp[,c("5084-0.0","5085-0.0","5201-0.0")]
UKBB.phenotypes.temp1[UKBB.phenotypes.temp1<0] <- NA
UKBB.phenotypes.temp <- cbind(UKBB.phenotypes.temp1,UKBB.phenotypes.temp2)
UKBB.phenotypes.temp[,benefit:=fifelse(`6146-0.0`==100,0,1,na=NA)]
UKBB.phenotypes.temp[,degree_college:=fifelse(`6138-0.0`==1,1,0,na=NA)]
UKBB.phenotypes.temp[,degree_degree:=fifelse(`6138-0.0`==100,0,1,na=NA)]
UKBB.phenotypes.temp[`4825-0.0`==0,`4825-0.0`:=14]
UKBB.phenotypes.temp[`1100-0.0`==5,`1100-0.0`:=NA]
UKBB.phenotypes.temp <- UKBB.phenotypes.temp[,-c("6146-0.0","6138-0.0")]
message(" ... done")
UKBB.phenotypes.parsed <- UKBB.phenotypes.temp

# Add FIS measures to cleaned imputation phenotype data
UKBB.phenotypes.parsed <- merge(UKBB.phenotypes.parsed, FIS_measures, by="eid", all.y=TRUE)

# Remove obsolete data 
rm(UKBB.phenotypes.temp)
rm(UKBB.phenotypes.temp1) 
rm(UKBB.phenotypes.temp2)
rm(UKBB.phenotypes)
rm(FIS_measures)

save.image(file = paste0(output_folder, "imputation_workspace_2.RData"))

# ---- 3. PERFORM IMPUTATION ---- #

# Keep original UKBB.phenotypes for reference
UKBB.phenotypes <- UKBB.phenotypes.parsed

message("3. Performing imputation")
UKBB.phenotypes.mat <- as.matrix(UKBB.phenotypes)
UKBB.phenotypes.mat <- matrix(as.numeric(UKBB.phenotypes.mat),ncol=ncol(UKBB.phenotypes.mat))
UKBB.phenotypes.mat.s <- biScale(UKBB.phenotypes.mat)
# Rank and lambda are 150, 120 for v1. rest is 80, 70
fit=softImpute(UKBB.phenotypes.mat.s,rank=80,lambda=70) 
UKBB.phenotypes.imputed=complete(UKBB.phenotypes.mat,fit)
UKBB.phenotypes.imputed <- as.data.table(UKBB.phenotypes.imputed)

setnames(UKBB.phenotypes.imputed,names(UKBB.phenotypes))
message(" ... done")

rm(UKBB.phenotypes)
rm(UKBB.phenotypes.mat)
rm(UKBB.phenotypes.mat.s)

save.image(file = paste0(output_folder, "imputation_workspace_3.RData"))

# ---- 4. PREPARE PHENO FOR GWAS ---- #
message("4. Preparing output")
# Subset measured and imputed FIS values, change (20016-0.0, mean_FIS, earliest) accordingly
UKBB.FIS.all <- UKBB.phenotypes.imputed[,c('eid', 'eid', '20016-0.0')]
colnames(UKBB.FIS.all) <- c("EID", "FID", "PHENO")

UKBB.FIS.measured <- UKBB.phenotypes.parsed[,c('eid', 'eid', '20016-0.0')]
UKBB.FIS.measured <- na.omit(UKBB.FIS.measured)
colnames(UKBB.FIS.measured) <- c("EID", "FID", "PHENO")

UKBB.FIS.imputed <- anti_join(UKBB.FIS.all, UKBB.FIS.measured, by = "EID")
colnames(UKBB.FIS.imputed) <- c("EID", "FID", "PHENO")

#Function to remove extreme outliers from imputed FIS values
filter_outliers <- function(imp_df, obs_df, col_name) {
  obs_mean_val <- mean(obs_df[[col_name]], na.rm = TRUE)
  obs_sd_val <- sd(obs_df[[col_name]], na.rm = TRUE)
  df_outliers <- imp_df[abs(imp_df[[col_name]] - obs_mean_val) <= 3 * obs_sd_val, ]
  return(df_outliers)
}

initial_imputed <- UKBB.FIS.imputed$"EID"
UKBB.FIS.imputed <- filter_outliers(UKBB.FIS.imputed, UKBB.FIS.measured, 'PHENO')
message(paste0(" ... removed ", length(initial_imputed)-length(UKBB.FIS.imputed$'EID'), " extreme outliers"))
message(paste(setdiff(initial_imputed, UKBB.FIS.imputed$'EID'), collapse = "\n- "))
#Also remove outliers from all set
UKBB.FIS.all <- UKBB.FIS.all[!UKBB.FIS.all$EID %in% setdiff(initial_imputed, UKBB.FIS.imputed$'EID'), ]

message("\n --- OVERVIEW --- \n")
message(paste0(" ", length(colnames(UKBB.phenotypes.imputed))-1, " phenotypes used for imputation (including selected FIS measure(s))"))
message(paste0(" N measured: ", length(UKBB.FIS.measured$EID)))
message(paste0(" N imputed: ", length(UKBB.FIS.imputed$EID)))
message(paste0(" N all: ", length(UKBB.FIS.all$EID)))
save.image(file = paste0(output_folder, "imputation_workspace_4.RData"))

# Write files to pheno file
write.table(UKBB.FIS.imputed, paste0(output_folder,"imputed_FIS.pheno"), row.names = FALSE, sep='\t')
write.table(UKBB.FIS.measured, paste0(output_folder,"measured_FIS.pheno"), row.names = FALSE, sep='\t')

# --- Bonus: Create imputed ever measured and imputed never measured files --- #
never_measured <- UKBB.phenotypes.parsed %>% filter(across(all_of(included_FIS), is.na))
colnames(never_measured)[1] <- 'EID'
UKBB.FIS.imputed.ever <- anti_join(UKBB.FIS.imputed, never_measured, by="EID")
UKBB.FIS.imputed.never <- anti_join(UKBB.FIS.imputed, UKBB.FIS.imputed.ever, by="EID")
write.table(UKBB.FIS.imputed.ever, paste0(output_folder,"imputed_FIS_ever.pheno"), row.names = FALSE, sep='\t')
write.table(UKBB.FIS.imputed.never, paste0(output_folder,"imputed_FIS_never.pheno"), row.names = FALSE, sep='\t')
# ---------------------------------------------------------------------------- #
