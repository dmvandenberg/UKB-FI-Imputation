# GenmicSEM Cholesky model for IQ and ImputedIQ GWAS
# UPDATED VERSION for new imputation

require(GenomicSEM)
require(gdata)
require(Matrix)
require(parallel)

# loads LDSCoutput and sumstats
names = c("CTG_FIS","IMP_FIS")
load('LDSCoutputFISimpFIS.RData')
load('Sumstats.RData')


###### Cholesky model with SNP effect of F1 ##############################

#toevoeging Dirk Smit
source("userGWAS_dirk.R")
# alle userGWAS vervangen door userGWAS_dirk

# necessary functions for parallel processing. AddSNPs has been removed
fx_userGWAS <- function(a_sumstats_chunk, an_ldsc, 
                        modelGWAS, ncores=8, 
                        prefix='res', postfix="_BIPOCD") { 
  #stond eerst ncores=1
  ptm <- proc.time()
  FN = sprintf("%s%d%s.RData", prefix, a_sumstats_chunk$id[1], postfix)
  if (!file.exists(FN)) { #uitroepteken = nee / False
    cat(sprintf("Start %d\n",a_sumstats_chunk$id[1]))
    if (ncores>1) { #hiermee geeft ie aan of je kunt parallel processen
      #y <- userGWAS(covstruc=an_ldsc, SNPs=a_sumstats_chunk, model=modelGWAS, estimation = "DWLS", 
      #cores=ncores, parallel=T)
      y <- userGWAS(covstruc=an_ldsc,SNPs=a_sumstats_chunk, 
                         model=modelGWAS, estimation = "DWLS", 
                         cores=ncores, parallel=T)
    } else {
      #y <- userGWAS(covstruc=an_ldsc, SNPs=a_sumstats_chunk, model=modelGWAS, estimation = "DWLS", 
      #cores=1, parallel=F)
      y <- userGWAS(covstruc=an_ldsc, SNPs=a_sumstats_chunk, 
                         model=modelGWAS, estimation = "DWLS", 
                         cores=1, parallel=F)
    }
    save(file=FN, list = c("y"))
    cat(sprintf("Done %d\n",a_sumstats_chunk$id[1]))#wat geeft ie aan als ie done is
  } else {
    cat(sprintf("Skip %d\n",a_sumstats_chunk$id[1])) #als ie geskipt wordt
  }
  return (proc.time()-ptm)[3] # returns elapsed time
}



#------------------------------------------------------------------

steps=500 #in chunks
len=dim(p_sumstats)[1]
starts = round(seq(1,len+1,by=(len/(steps))))
stops = starts[2:(steps+1)]-1
starts = starts[1:steps]
p_sumstats_chunks = list()
for (i in 1:steps) {
  p_sumstats_chunks[[i]] = p_sumstats[starts[i]:stops[i],]
  p_sumstats_chunks[[i]]$id = i
}

#-------------------------------------------------------------------

# the model contains two factors (unique and shared)

ModelCholGWAS <- "FIS=~NA*IMP_FIS+CTG_FIS
NFIS=~NA*IMP_FIS
FIS~SNP
NFIS~SNP
FIS~~1*FIS
NFIS~~1*NFIS
FIS~~0*NFIS
IMP_FIS~~0*IMP_FIS
CTG_FIS~~0*CTG_FIS
CTG_FIS~~0*IMP_FIS
SNP~~SNP"

# to test, run with one or two chunks: p_sumstats_chunks[1:2]
ptm <- proc.time()
par_results_userGWAS <- mclapply(X = p_sumstats_chunks, FUN = fx_userGWAS, 
                               an_ldsc=LDSCoutput, modelGWAS=ModelCholGWAS, 
                               ncores=1, 
                               prefix='res', postfix='_FIS_NFIS', mc.cores=8) 
# writes models for rsIDs to harddrive
proc.time()-ptm




