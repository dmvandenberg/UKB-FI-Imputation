#userGWAS_dirk <- function (Output, estimation = "DWLS", model = "", modelchi = FALSE, #origineel
userGWAS_dirk<- function (covstruc = NULL, SNPs = NULL, estimation = "DWLS", 
                          model = "", printwarn = TRUE, sub = FALSE, cores = NULL, 
                          toler = FALSE, SNPSE = FALSE, parallel = TRUE, GC = "standard", 
                          MPI = FALSE, smooth_check = FALSE, TWAS = FALSE, std.lv = FALSE, 
                          fix_measurement = TRUE) 
{
  require(GenomicSEM)
  require(lavaan)  
  require(stringr)
  require(splitstackshape)
  require(gdata)
  require(usethis) 
  require(doParallel)
  
  
  if (!toler) 
    toler <- .Machine$double.eps
  GenomicSEM:::.check_one_of(estimation, c("DWLS", "ML"))
  GenomicSEM:::.check_boolean(printwarn)
  if (!is.null(cores)) 
    GenomicSEM:::.check_range(cores, min = 0, max = Inf)
  GenomicSEM:::.check_range(toler, min = 0, max = Inf)
  GenomicSEM:::.check_boolean(parallel)
  GenomicSEM:::.check_one_of(GC, c("standard", "conserv", "none"))
  GenomicSEM:::.check_boolean(MPI)
  GenomicSEM:::.check_boolean(smooth_check)
  GenomicSEM:::.check_boolean(TWAS)
  GenomicSEM:::.check_boolean(std.lv)
  time <- proc.time()
  Operating <- Sys.info()[["sysname"]]
  if (MPI == TRUE & Operating == "Windows") {
    stop("MPI is not currently available for Windows operating systems. Please set the MPI argument to FALSE, or switch to a Linux or Mac operating system.")
  }
  if (exists("Output")) {
    stop("Please note that an update was made to commonfactorGWAS on 4/1/21 so that addSNPs output CANNOT be fed directly to the function. It now expects the\n            output from ldsc (using covstruc = ...)  followed by the output from sumstats (using SNPs = ... ) as the first two arguments.")
  }
  test <- c(str_detect(model, "~"), str_detect(model, "="), 
            str_detect(model, "\\+"))
  if (!all(test)) {
    warning("Your model name may be listed in quotes; please remove the quotes and try re-running if the function has returned stopped running after returning an error.")
  }
  print("Please note that an update was made to userGWAS on 11/21/19 so that it combines addSNPs and userGWAS.")
  if (class(SNPs)[1] == "character") {
    print("You are likely listing arguments in the order of a previous version of userGWAS, if you have yur results stored after running addSNPs you can still explicitly call Output = ... to provide them to userGWAS. The current version of the function is faster and saves memory. It expects the\n          output from ldsc followed by the output from sumstats (using SNPs = ... ) as the first two arguments. See ?userGWAS for help on propper usag")
    warning("You are likely listing arguments (e.g. Output = ...) in the order of a previous version of userGWAS, if you have yur results stored after running addSNPs you can still explicitly call Output = ... to provide them to userGWAS. The current version of the function is faster and saves memory. It expects the\n            output from ldsc (using covstruc = ...)  followed by the output from sumstats (using SNPs = ... ) as the first two arguments. See ?userGWAS for help on propper usage")
  }
  else {
    if (is.null(SNPs) | is.null(covstruc)) {
      print("You may be listing arguments in the order of a previous version of userGWAS, if you have yur results stored after running addSNPs you can still explicitly call Output = ... to provide them to userGWAS;if you already did this and the function ran then you can disregard this warning. The current version of the function is faster and saves memory. It expects the\n            output from ldsc followed by the output from sumstats (using SNPs = ... ) as the first two arguments. See ?userGWAS for help on propper usag")
      warning("You may be listing arguments (e.g. Output = ...) in the order of a previous version of userGWAS, if you have yur results stored after running addSNPs you can still explicitly call Output = ... to provide them to userGWAS; ; if you already did this and the function ran then you can disregard this warning. The current version of the function is faster and saves memory. It expects the\n              output from ldsc (using covstruc = ...)  followed by the output from sumstats (using SNPs = ... ) as the first two arguments. See ?userGWAS for help on propper usage")
    }
  }
  if (sub[[1]] != FALSE) {
    sub <- str_replace_all(sub, fixed(" "), "")
  }
  SNPs <- data.frame(SNPs)
  if (TWAS) {
    SNPs$Gene <- as.character(SNPs$Gene)
    SNPs$Panel <- as.character(SNPs$Panel)
    varSNP <- SNPs$HSQ
  }
  else {
    SNPs$A1 <- as.character(SNPs$A1)
    SNPs$A2 <- as.character(SNPs$A2)
    SNPs$SNP <- as.character(SNPs$SNP)
    varSNP <- 2 * SNPs$MAF * (1 - SNPs$MAF)
  }
  if (SNPSE == FALSE) {
    varSNPSE2 <- (5e-04)^2
  }
  if (SNPSE != FALSE) {
    varSNPSE2 <- SNPSE^2
  }
  V_LD <- as.matrix(covstruc[[1]])
  S_LD <- as.matrix(covstruc[[2]])
  I_LD <- as.matrix(covstruc[[3]])
  Model1 <- model
  for (i in 1) {
    if (fix_measurement) {
      rownames(S_LD) <- colnames(S_LD)
      lines <- strsplit(model, "\n")[[1]]
      filtered_lines <- lines[!grepl("SNP", lines)]
      noSNPmodel <- paste(filtered_lines, collapse = "\n")
      smoothS <- ifelse(eigen(S_LD)$values[nrow(S_LD)] <= 
                          0, S_LD <- as.matrix((nearPD(S_LD, corr = FALSE))$mat), 
                        S_LD <- S_LD)
      smoothV <- ifelse(eigen(V_LD)$values[nrow(V_LD)] <= 
                          0, V_LD <- as.matrix((nearPD(V_LD, corr = FALSE))$mat), 
                        V_LD <- V_LD)
      W <- solve(V_LD, tol = toler)
      testnoSNP <- GenomicSEM:::.tryCatch.W.E(ReorderModelnoSNP <- sem(noSNPmodel, 
                                                          sample.cov = S_LD, estimator = "DWLS", WLS.V = W, 
                                                          sample.nobs = 2, optim.dx.tol = +Inf, optim.force.converged = TRUE, 
                                                          control = list(iter.max = 1), std.lv = std.lv))
      order <- GenomicSEM:::.rearrange(k = ncol(S_LD), fit = ReorderModelnoSNP, 
                          names = colnames(S_LD))
      V_Reorder <- V_LD[order, order]
      u <- nrow(V_Reorder)
      W_Reorder <- diag(u)
      diag(W_Reorder) <- diag(V_Reorder)
      W_Reorder <- solve(W_Reorder, tol = toler)
      if (estimation == "DWLS") {
        if (std.lv == FALSE) {
          emptynoSNP <- GenomicSEM:::.tryCatch.W.E(Model1_Results <- sem(noSNPmodel, 
                                                            sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, 
                                                            sample.nobs = 2, optim.dx.tol = +Inf))
        }
        if (std.lv == TRUE) {
          emptynoSNP <- GenomicSEM:::.tryCatch.W.E(Model1_Results <- sem(noSNPmodel, 
                                                            sample.cov = S_LD, estimator = "DWLS", WLS.V = W_Reorder, 
                                                            sample.nobs = 2, std.lv = TRUE, optim.dx.tol = +Inf))
        }
      }
      if (estimation == "ML") {
        if (std.lv == FALSE) {
          emptynoSNP <- GenomicSEM:::.tryCatch.W.E(Model1_Results <- sem(noSNPmodel, 
                                                            sample.cov = S_LD, estimator = "ML", sample.nobs = 200, 
                                                            optim.dx.tol = +Inf, sample.cov.rescale = FALSE))
        }
        if (std.lv == TRUE) {
          emptynoSNP <- GenomicSEM:::.tryCatch.W.E(Model1_Results <- sem(noSNPmodel, 
                                                            sample.cov = S_LD, estimator = "ML", sample.nobs = 200, 
                                                            std.lv = TRUE, optim.dx.tol = +Inf, sample.cov.rescale = FALSE))
        }
      }
      Model1 <- parTable(Model1_Results)
      for (p in 1:nrow(Model1)) {
        Model1$free[p] <- ifelse(Model1$lhs[p] != Model1$rhs[p], 
                                 0, Model1$free[p])
      }
    }
  }
  beta_SNP <- SNPs[, grep("beta.", fixed = TRUE, colnames(SNPs))]
  SE_SNP <- SNPs[, grep("se.", fixed = TRUE, colnames(SNPs))]
  n_phenotypes <- ncol(beta_SNP)
  diag(I_LD) <- ifelse(diag(I_LD) <= 1, 1, diag(I_LD))
  coords <- which(I_LD != "NA", arr.ind = T)
  i <- 1
  V_SNP <- diag(n_phenotypes)
  for (p in 1:nrow(coords)) {
    x <- coords[p, 1]
    y <- coords[p, 2]
    if (x != y) {
      V_SNP[x, y] <- (SE_SNP[i, y] * SE_SNP[i, x] * I_LD[x, 
                                                         y] * I_LD[x, x] * I_LD[y, y] * varSNP[i]^2)
    }
    if (x == y) {
      V_SNP[x, x] <- (SE_SNP[i, x] * I_LD[x, x] * varSNP[i])^2
    }
  }
  V_full <- GenomicSEM:::.get_V_full(n_phenotypes, V_LD, varSNPSE2, V_SNP)
  kv <- nrow(V_full)
  smooth2 <- ifelse(eigen(V_full)$values[kv] <= 0, V_full <- as.matrix((nearPD(V_full, 
                                                                               corr = FALSE))$mat), V_full <- V_full)
  S_SNP <- vector(mode = "numeric", length = n_phenotypes + 
                    1)
  S_SNP[1] <- varSNP[i]
  for (p in 1:n_phenotypes) {
    S_SNP[p + 1] <- varSNP[i] * beta_SNP[i, p]
  }
  S_Full <- diag(n_phenotypes + 1)
  S_Full[(2:(n_phenotypes + 1)), (2:(n_phenotypes + 1))] <- S_LD
  S_Full[1:(n_phenotypes + 1), 1] <- S_SNP
  S_Full[1, 1:(n_phenotypes + 1)] <- t(S_SNP)
  if (TWAS) {
    colnames(S_Full) <- c("Gene", colnames(S_LD))
  }
  else {
    colnames(S_Full) <- c("SNP", colnames(S_LD))
  }
  rownames(S_Full) <- colnames(S_Full)
  ks <- nrow(S_Full)
  smooth1 <- ifelse(eigen(S_Full)$values[ks] <= 0, S_Full <- as.matrix((nearPD(S_Full, 
                                                                               corr = FALSE))$mat), S_Full <- S_Full)
  k2 <- ncol(S_Full)
  for (i in 1) {
    W <- solve(V_full, tol = toler)
    test2 <- GenomicSEM:::.tryCatch.W.E(ReorderModel <- sem(model, sample.cov = S_Full, 
                                               estimator = "DWLS", WLS.V = W, sample.nobs = 2, optim.dx.tol = +Inf, 
                                               optim.force.converged = TRUE, control = list(iter.max = 1), 
                                               std.lv = std.lv))
    if (fix_measurement) {
      withSNP <- parTable(ReorderModel)
      for (p in 1:nrow(withSNP)) {
        if (withSNP$rhs[p] == "SNP" | withSNP$lhs[p] == 
            "SNP") {
          Model1 <- rbind(Model1, withSNP[p, ])
        }
      }
      test3 <- GenomicSEM:::.tryCatch.W.E(ReorderModel <- sem(Model1, 
                                                 sample.cov = S_Full, estimator = "DWLS", WLS.V = W, 
                                                 sample.nobs = 2, optim.dx.tol = +Inf, optim.force.converged = TRUE, 
                                                 control = list(iter.max = 1), std.lv = std.lv))
    }
    order <- GenomicSEM:::.rearrange(k = k2, fit = ReorderModel, names = rownames(S_Full))
    suppressWarnings(df <- lavInspect(ReorderModel, "fit")["df"])
    suppressWarnings(npar <- lavInspect(ReorderModel, "fit")["npar"])
  }
  if (TWAS) {
    SNPs2 <- SNPs[, 1:3]
  }
  else {
    SNPs2 <- SNPs[, 1:6]
  }
  rm(SNPs)
  f <- nrow(beta_SNP)
  LavModel1 <- GenomicSEM:::.userGWAS_main(i = 1, cores = 1, n_phenotypes, 
                              n = 1, I_LD, V_LD, S_LD, std.lv, varSNPSE2, order, SNPs2, 
                              beta_SNP, SE_SNP, varSNP, GC, coords, smooth_check, TWAS, 
                              printwarn, toler, estimation, sub, Model1, df, npar, 
                              returnlavmodel = TRUE)
  if (!parallel) {
    if (sub[[1]] == FALSE) {
      Results_List <- vector(mode = "list", length = nrow(beta_SNP))
    }
    if (TWAS) {
      print("Starting TWAS Estimation")
    }
    else {
      print("Starting GWAS Estimation")
    }
    for (i in 1:nrow(beta_SNP)) {
      if (i == 1) {
        cat(paste0("Running Model: ", i, "\n"))
      }
      else {
        if (i%%1000 == 0) {
          cat(paste0("Running Model: ", i, "\n"))
        }
      }
      final2 <- GenomicSEM:::.userGWAS_main(i, cores = 1, n_phenotypes, 
                               1, I_LD, V_LD, S_LD, std.lv, varSNPSE2, order, 
                               SNPs2, beta_SNP, SE_SNP, varSNP, GC, coords, 
                               smooth_check, TWAS, printwarn, toler, estimation, 
                               sub, Model1, df, npar, basemodel = LavModel1)
      final2$i <- NULL
      if (sub[[1]] != FALSE) {
        final3 <- as.data.frame(matrix(NA, ncol = ncol(final2), 
                                       nrow = length(sub)))
        final3[1:length(sub), ] <- final2[1:length(sub), 
        ]
        if (i == 1) {
          Results_List <- vector(mode = "list", length = length(sub))
          for (y in 1:length(sub)) {
            Results_List[[y]] <- as.data.frame(matrix(NA, 
                                                      ncol = ncol(final3), nrow = f))
            colnames(Results_List[[y]]) <- colnames(final2)
            Results_List[[y]][1, ] <- final3[y, ]
          }
        }
        else {
          for (y in 1:nrow(final3)) {
            Results_List[[y]][i, ] <- final3[y, ]
          }
        }
      }
      else {
        Results_List[[i]] <- final2
      }
    }
    time_all <- proc.time() - time
    print(time_all[3])
    return(Results_List)
  }
  else {
    if (is.null(cores)) {
      int <- min(c(nrow(SNPs2), detectCores() - 1))
    }
    else {
      if (cores > nrow(SNPs2)) 
        warning(paste0("Provided number of cores was greater than number of SNPs, reverting to cores=", 
                       nrow(SNPs2)))
      int <- min(c(cores, nrow(SNPs2)))
    }
    if (MPI) {
      cl <- getMPIcluster()
      registerDoParallel(cl)
    }
    else {
      if (Operating != "Windows") {
        cl <- makeCluster(int, type = "FORK")
      }
      else {
        cl <- makeCluster(int, type = "PSOCK")
      }
      registerDoParallel(cl)
      on.exit(stopCluster(cl))
    }
    SNPs2 <- suppressWarnings(split(SNPs2, 1:int))
    beta_SNP <- suppressWarnings(split(beta_SNP, 1:int))
    SE_SNP <- suppressWarnings(split(SE_SNP, 1:int))
    varSNP <- suppressWarnings(split(varSNP, 1:int))
    if (TWAS) {
      print("Starting TWAS Estimation")
    }
    else {
      print("Starting GWAS Estimation")
    }
    if (Operating != "Windows") {
      results <- foreach(n = icount(int), .combine = "rbind") %:% 
        foreach(i = 1:nrow(beta_SNP[[n]]), .combine = "rbind", 
                .packages = "lavaan") %dopar% GenomicSEM:::.userGWAS_main(i, 
                                                             int, n_phenotypes, n, I_LD, V_LD, S_LD, std.lv, 
                                                             varSNPSE2, order, SNPs2[[n]], beta_SNP[[n]], 
                                                             SE_SNP[[n]], varSNP[[n]], GC, coords, smooth_check, 
                                                             TWAS, printwarn, toler, estimation, sub, Model1, 
                                                             df, npar, basemodel = LavModel1)
    }
    else {
      utilfuncs <- list()
      utilfuncs[["GenomicSEM:::.tryCatch.W.E"]] <- GenomicSEM:::.tryCatch.W.E
      utilfuncs[["GenomicSEM:::.get_V_SNP"]] <- GenomicSEM:::.get_V_SNP
      utilfuncs[["GenomicSEM:::.get_Z_pre"]] <- GenomicSEM:::.get_Z_pre
      utilfuncs[["GenomicSEM:::.get_V_full"]] <- GenomicSEM:::.get_V_full
      results <- foreach(n = icount(int), .combine = "rbind") %:% 
        foreach(i = 1:nrow(beta_SNP[[n]]), .combine = "rbind", 
                .packages = c("lavaan", "gdata"), .export = c(".userGWAS_main")) %dopar% 
        {
          GenomicSEM:::.userGWAS_main(i, int, n_phenotypes, n, I_LD, 
                         V_LD, S_LD, std.lv, varSNPSE2, order, SNPs2[[n]], 
                         beta_SNP[[n]], SE_SNP[[n]], varSNP[[n]], 
                         GC, coords, smooth_check, TWAS, printwarn, 
                         toler, estimation, sub, Model1, df, npar, 
                         utilfuncs, basemodel = LavModel1)
        }
    }
    results <- results[order(results$i), ]
    results$i <- NULL
    if (sub[[1]] != FALSE) {
      Results_List <- vector(mode = "list", length = length(sub))
      for (y in 1:length(sub)) {
        Results_List[[y]] <- as.data.frame(matrix(NA, 
                                                  ncol = ncol(results), nrow = nrow(results)/length(sub)))
        colnames(Results_List[[y]]) <- colnames(results)
        Results_List[[y]] <- subset(results, paste0(results$lhs, 
                                                    results$op, results$rhs, sep = "") %in% sub[[y]] | 
                                      is.na(results$lhs))
      }
      rm(results)
    }
    if (sub[[1]] == FALSE) {
      if (TWAS) {
        names <- unique(results$Panel)
        Results_List <- vector(mode = "list", length = length(names))
        for (y in 1:length(names)) {
          Results_List[[y]] <- subset(results, results$Panel == 
                                        names[[y]])
          Results_List[[y]]$Model_Number <- NULL
        }
      }
      else {
        names <- unique(results$SNP)
        Results_List <- vector(mode = "list", length = length(names))
        for (y in 1:length(names)) {
          Results_List[[y]] <- subset(results, results$SNP == 
                                        names[[y]])
          Results_List[[y]]$Model_Number <- NULL
        }
      }
      rm(results)
      rm(names)
    }
    time_all <- proc.time() - time
    print(time_all[3])
    return(Results_List)
  }
}
