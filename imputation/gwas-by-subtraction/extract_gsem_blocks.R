# Parse gwas-by-sub chunks
library(dplyr)

process_file <- function(file_path) {
  load(file_path)
  # Extract NFIS ~ SNP line
  nfis_snp <- lapply(y, function(x) x[11, ])
  # Combine the rows into a single data frame
  nfis_chunk <- do.call(rbind, nfis_snp)
  rm(y)
  print(paste0("Extraction complete: ", file_path))
  return(nfis_chunk)
}

gwas_res <- "./GWAS-by-sub"
rdata_files <- list.files(gwas_res, pattern = "\\.RData$", full.names = TRUE)
all_chunks <- lapply(rdata_files, process_file)
nfis_gwas <- bind_rows(all_chunks)

# Perform some checks
# 1. Correct lines extracted
unique(nfis_gwas$lhs)
unique(nfis_gwas$rhs)

# 2. No errors or warnings
unique(nfis_gwas$error)
unique(nfis_gwas$warning)

# Extract columns and rename
nfis_gwas_clean <- nfis_gwas[, c("SNP", "CHR", "BP", "A1", "A2", "MAF", "est", "SE", "Z_Estimate", "Pval_Estimate")]
colnames(nfis_gwas_clean) <- c("SNP", "CHR", "BP", "A1", "A2", "MAF", "BETA", "SE", "Z", "P")

# Remove SNPs with MAF <= .1 and MAF >= .4 (as recommended)
nfis_gwas_clean <- nfis_gwas_clean[nfis_gwas_clean$'MAF' >= .1 & nfis_gwas_clean$'MAF' <= .4,]
# Compute N_hat (as recommended)
N_hat<-mean(1/((2*nfis_gwas_clean$MAF*(1-nfis_gwas_clean$MAF))*nfis_gwas_clean$SE^2))
nfis_gwas_clean$N <- round(N_hat, 0)

# Write file
write.table(nfis_gwas_clean, "NonCogIFIS.pheno", sep="\t", row.names = FALSE, quote = FALSE)
