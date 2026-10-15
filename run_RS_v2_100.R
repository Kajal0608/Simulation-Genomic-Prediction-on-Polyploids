### =========================================================
###  Recurrent selection 
### =========================================================

suppressPackageStartupMessages({
  library(AlphaSimR)
  library(rrBLUP)
  library(data.table)
})

t_master <- Sys.time()
cat("Starting at:", format(t_master), "\n")
cat("AlphaSimR:", as.character(packageVersion("AlphaSimR")), "\n")

## =========================================================
## Functions
## =========================================================
call_Bmid <- function(M3){
  M2 <- M3
  M2[M3 == 2] <- 1
  M2[M3 == 3] <- 2
  storage.mode(M2) <- "integer"
  M2
}

project_KA <- function(Mtr, Mte = NULL, m){
  p <- colMeans(Mtr) / m
  Wtr <- scale(Mtr, center = m * p, scale = FALSE)
  denom <- m * sum(p * (1 - p))
  Ktr <- tcrossprod(Wtr) / denom
  if (is.null(Mte)) return(list(Ktr = as.matrix(Ktr)))
  Wte <- scale(Mte, center = m * p, scale = FALSE)
  Kte <- (Wte %*% t(Wtr)) / denom
  list(Ktr = as.matrix(Ktr), Kte = as.matrix(Kte))
}

fit_predict_KGBLUP <- function(ytr, Ktr, Kte){
  ybar <- mean(ytr)
  yctr <- ytr - ybar
  fit  <- rrBLUP::mixed.solve(y = ytr, K = Ktr)
  lam  <- as.numeric(fit$Ve / fit$Vu)
  a    <- solve(Ktr + diag(lam, nrow(Ktr)), yctr)
  yhat <- as.numeric(ybar + Kte %*% a)
  list(yhat = yhat, lam = lam)
}

## =========================================================
## 1) Pool initialization 
## =========================================================
init_pools <- function(
    seed_main = 29121983,
    nChr = 11, segSites = 3500,
    nQTL = 495, chip.size.1 = 3003,
    use_gamma = TRUE, varA_target = 1,
    meanDD = 0.6, varDD = 0.5, h2_base = 0.12,
    n_founders_hist = 1000,
    N2 = 400, N4 = 400, quadProb = 0.3
){
  set.seed(seed_main)
  qtl_per_chr <- nQTL / nChr
  snp_per_chr <- chip.size.1 / nChr
  segSites_per_chr <- as.integer(segSites / nChr)
  Ve_fixed <- varA_target * (1 - h2_base) / h2_base
  
  history2 <- quickHaplo(nInd = n_founders_hist, nChr = nChr, 
                          segSites = segSites_per_chr)
  
  SP2 <- SimParam$new(history2)
  SP2$restrSegSites(qtl_per_chr, snp_per_chr)
  SP2$addTraitAD(nQtlPerChr = qtl_per_chr,
                 mean = 0, var = varA_target, useVarA = TRUE,
                 meanDD = meanDD, varDD = varDD, gamma = use_gamma)
  SP2$setVarE(varE = Ve_fixed)
  SP2$addSnpChip(snp_per_chr)
  
  pop2 <- newPop(rawPop = history2, isDH = FALSE, simParam = SP2)
  pop2@id <- paste0("D0_", sprintf("%05d", 1:nInd(pop2)))
  pool2 <- pop2[sample(1:nInd(pop2), N2)]
  
  base2_for4 <- pop2[sample(1:nInd(pop2), max(N4, 250))]
  pop4_tmp <- doubleGenome(base2_for4, simParam = SP2)
  
  H4 <- pullSegSiteHaplo(pop4_tmp, simParam = SP2)
  cols_per_chr <- ncol(H4) / nChr
  H4_list <- lapply(seq_len(nChr), function(i){
    a <- (i - 1) * cols_per_chr + 1
    b <- i * cols_per_chr
    H4[, a:b, drop = FALSE]
  })
  
  founders4 <- newMapPop(genMap = SP2$genMap, haplotypes = H4_list, ploidy = 4L)
  
  SP4 <- SimParam$new(founders4)
  SP4$restrSegSites(qtl_per_chr, snp_per_chr)
  SP4$manAddTrait(SP2$traits[[1]])
  SP4$setVarE(varE = Ve_fixed)
  SP4$addSnpChip(snp_per_chr)
  SP4$quadProb <- quadProb
  
  pop4 <- newPop(rawPop = founders4, simParam = SP4)
  pop4@id <- paste0("T0_", sprintf("%05d", 1:nInd(pop4)))
  pool4 <- pop4[sample(1:nInd(pop4), N4)]
  
  list(SP2 = SP2, SP4 = SP4, pool2 = pool2, pool4 = pool4, Ve_fixed = Ve_fixed)
}

## =========================================================
## 2) One RS cycle
## =========================================================
run_one_cycle <- function(
    cycle, SP2, SP4, pool2, pool4,
    calling = c("C", "Bmid"),
    n_f_train = 60, n_m_train = 60,
    n_cross_train = 300, n_prog_train = 5,
    k2 = 40, k4 = 40,
    N2_next = 400, N4_next = 400,
    n_f_eval = 40, n_m_eval = 40,
    n_cross_eval = 200, n_prog_eval = 3,
    seed_cross = 9001, seed_pheno = 111,
    train_target = c("TBV", "PHENO")
){
  calling <- match.arg(calling)
  train_target <- match.arg(train_target)
  
  ## ---- Training triploids ----
  set.seed(seed_cross + 1000 * cycle)
  fem_train <- pool4[sample(1:nInd(pool4), n_f_train)]
  mal_train <- pool2[sample(1:nInd(pool2), n_m_train)]
  
  cross_plan <- cbind(
    female = sample(1:n_f_train, n_cross_train, replace = TRUE),
    male   = sample(1:n_m_train, n_cross_train, replace = TRUE)
  )
  
  trPop <- makeCross2(fem_train, mal_train, crossPlan = cross_plan,
                      nProgeny = n_prog_train, simParam = SP4)
  
  set.seed(seed_pheno + 1000 * cycle)
  trPop <- setPheno(trPop, simParam = SP4)
  
  ytr <- if (train_target == "TBV") as.numeric(trPop@gv[, 1]) 
         else as.numeric(pheno(trPop)[, 1])
  M3_tr <- as.matrix(pullSnpGeno(trPop, simParam = SP4))
  
  Mtr   <- if (calling == "C") M3_tr else call_Bmid(M3_tr)
  m_use <- if (calling == "C") 3L else 2L
  
  P <- project_KA(Mtr, m = m_use)
  Ktr <- P$Ktr
  
 
  M4_parents <- as.matrix(pullSnpGeno(pool4, simParam = SP4))  # 0..4 dosage
  

  pool2_freq <- colMeans(as.matrix(pullSnpGeno(pool2, simParam = SP4))) / 2
  M3_exp_4 <- 0.5 * M4_parents + matrix(pool2_freq, 
                                         nrow = nInd(pool4),
                                         ncol = length(pool2_freq), 
                                         byrow = TRUE)
  
  Mpred_4 <- if (calling == "C") M3_exp_4 else pmin(round(M3_exp_4), 2)
  
  P4 <- project_KA(Mtr, Mpred_4, m = m_use)
  pred4 <- fit_predict_KGBLUP(ytr = ytr, Ktr = Ktr, Kte = P4$Kte)
  
  GCA4 <- data.table(parent_id = pool4@id, yhat = pred4$yhat)
  setorder(GCA4, -yhat)
  sel4_ids <- head(GCA4$parent_id, k4)
  

  M2_parents <- as.matrix(pullSnpGeno(pool2, simParam = SP4))  # 0..2 dosage
  

  pool4_freq <- colMeans(M4_parents) / 4
  M3_exp_2 <- M2_parents + matrix(2 * pool4_freq,
                                   nrow = nInd(pool2),
                                   ncol = length(pool4_freq),
                                   byrow = TRUE)
  
  Mpred_2 <- if (calling == "C") M3_exp_2 else pmin(round(M3_exp_2), 2)
  
  P2 <- project_KA(Mtr, Mpred_2, m = m_use)
  pred2 <- fit_predict_KGBLUP(ytr = ytr, Ktr = Ktr, Kte = P2$Kte)
  
  GCA2 <- data.table(parent_id = pool2@id, yhat = pred2$yhat)
  setorder(GCA2, -yhat)
  sel2_ids <- head(GCA2$parent_id, k2)
  

  sel4 <- pool4[match(sel4_ids, pool4@id)]
  sel2 <- pool2[match(sel2_ids, pool2@id)]
  
  next4 <- randCross(pop = sel4, nCrosses = ceiling(N4_next / 2), 
                     nProgeny = 2, simParam = SP4)
  next2 <- randCross(pop = sel2, nCrosses = ceiling(N2_next / 2),
                     nProgeny = 2, simParam = SP2)
  
  if (nInd(next4) > N4_next) next4 <- next4[1:N4_next]
  if (nInd(next2) > N2_next) next2 <- next2[1:N2_next]
  
  next4@id <- paste0("T", cycle, "_", sprintf("%05d", 1:nInd(next4)))
  next2@id <- paste0("D", cycle, "_", sprintf("%05d", 1:nInd(next2)))
  
  
  set.seed(seed_cross + 3000 * cycle)
  fem_eval <- next4[sample(1:nInd(next4), n_f_eval)]
  mal_eval <- next2[sample(1:nInd(next2), n_m_eval)]
  
  plan_eval <- cbind(
    female = sample(1:n_f_eval, n_cross_eval, replace = TRUE),
    male   = sample(1:n_m_eval, n_cross_eval, replace = TRUE)
  )
  
  evPop <- makeCross2(fem_eval, mal_eval, crossPlan = plan_eval,
                      nProgeny = n_prog_eval, simParam = SP4)
  
  mean_hyb_TBV <- mean(evPop@gv[, 1])
  mean_pool4 <- mean(next4@gv[, 1])
  mean_pool2 <- mean(next2@gv[, 1])
  
  list(
    next2 = next2,
    next4 = next4,
    stats = data.table(
      cycle = cycle, calling = calling,
      mean_hyb_TBV = mean_hyb_TBV,
      mean_pool4_TBV = mean_pool4,
      mean_pool2_TBV = mean_pool2
    )
  )
}

## =========================================================
## 3) Run one full RS 
## =========================================================
run_RS <- function(calling, T = 10, init_args = list(), cycle_args = list()){
  init <- do.call(init_pools, init_args)
  SP2 <- init$SP2; SP4 <- init$SP4
  pool2 <- init$pool2; pool4 <- init$pool4
  
  out_stats <- vector("list", T)
  for (t in 1:T) {
    t0_cycle <- Sys.time()
    step <- do.call(run_one_cycle, c(list(
      cycle = t, SP2 = SP2, SP4 = SP4, 
      pool2 = pool2, pool4 = pool4, calling = calling
    ), cycle_args))
    pool2 <- step$next2
    pool4 <- step$next4
    out_stats[[t]] <- step$stats
    cat(sprintf("    [%s] cycle %d/%d (%s) done in %.1f sec\n",
                format(Sys.time(), "%H:%M:%S"), t, T, calling,
                as.numeric(difftime(Sys.time(), t0_cycle, units = "secs"))))
  }
  rbindlist(out_stats)
}

## =========================================================
## 4) Replicate 
## =========================================================
run_RS_nreps <- function(nRep = 20, T = 10, init_args = list(), 
                         cycle_args = list(), seed_rep0 = 12345){
  all_stats <- vector("list", nRep)
  
  for (r in 1:nRep) {
    cat(sprintf("\n=== [%s] Rep %d / %d ===\n",
                format(Sys.time(), "%H:%M:%S"), r, nRep))
    
    init_args_r  <- modifyList(init_args, 
                               list(seed_main = seed_rep0 + 100000 * r))
    cycle_args_r <- modifyList(cycle_args, list(
      seed_cross = cycle_args$seed_cross + 10000 * r,
      seed_pheno = cycle_args$seed_pheno + 20000 * r
    ))
    
    st_C    <- run_RS("C",    T = T, init_args = init_args_r, 
                      cycle_args = cycle_args_r)
    st_Bmid <- run_RS("Bmid", T = T, init_args = init_args_r,
                      cycle_args = cycle_args_r)
    
    st <- rbind(st_C, st_Bmid)
    st[, rep := r]
    all_stats[[r]] <- st
  }
  
  out <- rbindlist(all_stats)
  setcolorder(out, c("rep", "cycle", "calling", "mean_hyb_TBV",
                     "mean_pool4_TBV", "mean_pool2_TBV"))
  out
}

## =========================================================
## 5) RUN
## =========================================================
init_args <- list(
  meanDD = 0.6, varDD = 0.5, h2_base = 0.12,
  N2 = 400, N4 = 400, quadProb = 0.3,
  n_founders_hist = 1000
)

cycle_args <- list(
  n_f_train = 60, n_m_train = 60,
  n_cross_train = 300, n_prog_train = 5,
  k2 = 40, k4 = 40, N2_next = 400, N4_next = 400,
  n_f_eval = 40, n_m_eval = 40,
  n_cross_eval = 200, n_prog_eval = 3,
  seed_cross = 9001, seed_pheno = 111,
  train_target = "TBV"
)

nRep_run <- 100
T_run    <- 10
seed_rep0 <- 777

cat(sprintf("\nLaunching RS: %d reps x %d cycles (periodic saves every 10 reps)\n",
            nRep_run, T_run))

all_stats <- vector("list", nRep_run)

for (r in 1:nRep_run) {
  cat(sprintf("\n=== [%s] Rep %d / %d ===\n",
              format(Sys.time(), "%H:%M:%S"), r, nRep_run))
  
  init_args_r  <- modifyList(init_args, 
                             list(seed_main = seed_rep0 + 100000 * r))
  cycle_args_r <- modifyList(cycle_args, list(
    seed_cross = cycle_args$seed_cross + 10000 * r,
    seed_pheno = cycle_args$seed_pheno + 20000 * r
  ))
  
  st_C    <- run_RS("C",    T = T_run, init_args = init_args_r, 
                    cycle_args = cycle_args_r)
  st_Bmid <- run_RS("Bmid", T = T_run, init_args = init_args_r,
                    cycle_args = cycle_args_r)
  
  st <- rbind(st_C, st_Bmid)
  st[, rep := r]
  all_stats[[r]] <- st
  
  if (r %% 10 == 0 || r == nRep_run) {
    partial <- rbindlist(all_stats[1:r])
    setcolorder(partial, c("rep", "cycle", "calling", "mean_hyb_TBV",
                           "mean_pool4_TBV", "mean_pool2_TBV"))
    saveRDS(partial, "results/stats_RS_v2_full100.rds")
    data.table::fwrite(partial, "results/stats_RS_v2_full100.csv")
    cat(sprintf("    >>> Saved partial: %d reps complete <<<\n", r))
  }
}

stats_all <- rbindlist(all_stats)
setcolorder(stats_all, c("rep", "cycle", "calling", "mean_hyb_TBV",
                         "mean_pool4_TBV", "mean_pool2_TBV"))

elapsed <- as.numeric(difftime(Sys.time(), t_master, units = "mins"))
cat(sprintf("\n=== RS DONE in %.1f min ===\n", elapsed))
cat(sprintf("Per cycle: %.2f sec\n", (elapsed * 60) / (nRep_run * T_run * 2)))

saveRDS(stats_all, "results/stats_RS_v2_full100.rds")
data.table::fwrite(stats_all, "results/stats_RS_v2_full100.csv")
cat("=== FINAL SAVED ===\n")

end_stats <- stats_all[cycle == T_run, 
                       .(mean_hyb_TBV = mean(mean_hyb_TBV)), 
                       by = .(rep, calling)]
end_wide <- dcast(end_stats, rep ~ calling, value.var = "mean_hyb_TBV")
end_wide[, diff_C_minus_Bmid := C - Bmid]

cat("\n=== Final-cycle summary across reps ===\n")
cat(sprintf("Mean C    final TBV: %.4f\n", mean(end_wide$C)))
cat(sprintf("Mean Bmid final TBV: %.4f\n", mean(end_wide$Bmid)))
cat(sprintf("Mean diff (C - Bmid): %.4f\n", mean(end_wide$diff_C_minus_Bmid)))
cat(sprintf("SE of diff: %.4f\n", 
            sd(end_wide$diff_C_minus_Bmid) / sqrt(nrow(end_wide))))

tt <- t.test(end_wide$C, end_wide$Bmid, paired = TRUE)
cat(sprintf("Paired t-test: t = %.3f, df = %d, p = %.4g\n",
            tt$statistic, tt$parameter, tt$p.value))
