### =========================================================
### Sensitivity study — when does dosage resolution matter? 

suppressPackageStartupMessages({
  library(AlphaSimR)
  library(rrBLUP)
  library(data.table)
})

## ----------------------------
## 0) Kernels + grouped CV
## ----------------------------
K_add <- function(M, m){
  pA <- colMeans(M) / m
  W  <- scale(M, center = m * pA, scale = FALSE)
  denom <- m * sum(pA * (1 - pA))
  tcrossprod(W) / denom
}

K_dom_trip <- function(M3){
  H <- 1.0 * ((M3 == 1) | (M3 == 2))
  q <- colMeans(H)
  Z <- scale(H, center=q, scale=FALSE)
  denom <- sum(q * (1 - q))
  tcrossprod(Z) / denom
}

K_dom_dip <- function(M2){
  stopifnot(all(M2 %in% 0:2))
  H <- 1.0 * (M2 == 1)
  q <- colMeans(H)
  Z <- scale(H, center=q, scale=FALSE)
  denom <- sum(q * (1 - q))
  tcrossprod(Z) / denom
}

cv_grouped_K <- function(y, Klist, groups, Kfold=5, seed=1){
  set.seed(seed)
  ug <- sample(unique(groups))
  folds <- split(ug, rep(1:Kfold, length.out=length(ug)))
  acc <- numeric(Kfold)
  
  for(k in seq_len(Kfold)){
    te_g <- folds[[k]]
    te <- which(groups %in% te_g)
    tr <- setdiff(seq_along(groups), te)
    
    ytr <- y[tr]; yte <- y[te]
    ybar <- mean(ytr)
    yctr <- ytr - ybar
    
    Ktr <- Reduce(`+`, lapply(Klist, function(Km) Km[tr, tr, drop=FALSE]))
    fit <- rrBLUP::mixed.solve(y=ytr, K=Ktr)
    lam <- as.numeric(fit$Ve / fit$Vu)
    
    a <- solve(Ktr + diag(lam, nrow(Ktr)), yctr)
    
    Kte <- Reduce(`+`, lapply(Klist, function(Km) Km[te, tr, drop=FALSE]))
    yhat <- as.numeric(ybar + Kte %*% a)
    
    acc[k] <- suppressWarnings(cor(yte, yhat, use="complete.obs"))
  }
  acc
}

## ----------------------------
## 1) One scenario runner 
## ----------------------------
run_scenario <- function(
    seed_main   = 29121983,
    seed_pheno  = 111,
    seed_cross  = 9001,
    seed_cv     = 123,
    
    # genome + architecture
    nChr        = 11,
    segSites    = 3500,
    nQTL        = 495,
    chip.size.1 = 3003,
    use_gamma   = TRUE,
    varA_target = 1,
    meanDD      = 0.6,
    varDD       = 0.5,
    
    # noise
    h2_base     = 0.12,
    
    # parents + selection
    n_founders_hist = 1000,
    n_parents_sel   = 96,
    
    # triploid training size knobs
    n_fem   = 40,
    n_mal   = 40,
    n_cross = 400,
    n_prog  = 5,
    
    # tetraploid pairing background (keep fixed for now)
    quadProb = 0.3,
    
    # CV
    Kfold  = 5
){
  
  set.seed(seed_main)
  
  ## derived counts
  qtl_per_chr <- nQTL / nChr
  snp_per_chr <- chip.size.1 / nChr
  stopifnot(qtl_per_chr == floor(qtl_per_chr),
            snp_per_chr == floor(snp_per_chr))
  segSites_per_chr <- as.integer(segSites / nChr)
  
  ## fixed Ve from target h2 in the *diploid founder scaling*
  Ve_fixed <- varA_target * (1 - h2_base) / h2_base
  
  ## ----- history haplotypes -----
  history2 <- quickHaplo(nInd = n_founders_hist, nChr = nChr, segSites = segSites_per_chr)
  
  ## ----- SP2 (diploid) with locked architecture -----
  SP2 <- SimParam$new(history2)
  SP2$restrSegSites(qtl_per_chr, snp_per_chr)
  SP2$addTraitAD(
    nQtlPerChr = qtl_per_chr,
    mean       = 0,
    var        = varA_target,
    meanDD     = meanDD,
    varDD      = varDD,
    useVarA    = TRUE,
    gamma      = use_gamma
  )
  SP2$setVarE(varE = Ve_fixed)
  SP2$addSnpChip(snp_per_chr)
  
  ## ----- diploid base pop + selection -----
  firstPop <- newPop(rawPop = history2, isDH = FALSE, simParam = SP2)
  set.seed(seed_pheno)
  firstPop <- setPheno(firstPop, simParam = SP2)
  
  firstPop_sel <- selectInd(firstPop,
                            nInd = n_parents_sel, trait = 1,
                            use = "pheno", gender = "B",
                            selectTop = TRUE,
                            returnPop = TRUE, simParam = SP2)
  
  ## ----- build tetraploids from elite diploids (doubleGenome) -----
  genMap2 <- tryCatch(SP2$genMap, error=function(e) NULL)
  if (is.null(genMap2)) genMap2 <- tryCatch(history2@genMap, error=function(e) NULL)
  stopifnot(!is.null(genMap2))
  
  firstPop4_tmp <- doubleGenome(firstPop_sel, simParam = SP2)
  H4 <- pullSegSiteHaplo(firstPop4_tmp, simParam = SP2)
  cols_per_chr <- ncol(H4) / nChr
  stopifnot(cols_per_chr == floor(cols_per_chr))
  
  H4_list <- lapply(seq_len(nChr), function(i){
    a <- (i-1)*cols_per_chr + 1
    b <- i*cols_per_chr
    H4[, a:b, drop=FALSE]
  })
  
  founders4_elite <- newMapPop(genMap = genMap2, haplotypes = H4_list, ploidy = 4L)
  
  SP4 <- SimParam$new(founders4_elite)
  SP4$restrSegSites(qtl_per_chr, snp_per_chr)
  SP4$manAddTrait(SP2$traits[[1]])
  SP4$setVarE(varE = Ve_fixed)
  SP4$addSnpChip(snp_per_chr)
  SP4$quadProb <- quadProb
  
  tetraPop_elite <- newPop(founders4_elite, simParam = SP4)
  
  ## ----- make triploid training population -----
  set.seed(seed_cross)
  
  femPop <- tetraPop_elite[sample(seq_len(nInd(tetraPop_elite)), n_fem)]
  malPop <- firstPop_sel[sample(seq_len(nInd(firstPop_sel)), n_mal)]
  
  cross_plan <- cbind(
    female = sample(1:n_fem, n_cross, replace=TRUE),
    male   = sample(1:n_mal, n_cross, replace=TRUE)
  )
  
  tripPop <- makeCross2(femPop, malPop, crossPlan = cross_plan, nProgeny = n_prog, simParam = SP4)
  tripPop <- setPheno(tripPop, simParam = SP4)
  
  ## dosages + TBV
  M3  <- as.matrix(pullSnpGeno(tripPop, simParam = SP4))  # 0..3
  tbv <- as.numeric(tripPop@gv[,1])
  
  cross_id <- rep(seq_len(nrow(cross_plan)), each = n_prog)
  stopifnot(length(cross_id) == nrow(M3))
  
  ## ----- C: proper triploid kernels -----
  KA_C <- K_add(M3, m=3)
  KD_C <- K_dom_trip(M3)
  
  ## ----- B_mid: low-cost diploid calling 0,1,2,3 -> 0,1,1,2 -----
  M2_mid <- M3
  M2_mid[M3 == 2] <- 1
  M2_mid[M3 == 3] <- 2
  stopifnot(all(M2_mid %in% 0:2))
  
  KA_B <- K_add(M2_mid, m=2)
  KD_B <- K_dom_dip(M2_mid)
  
  ## ----- grouped CV (by cross) -----
  set.seed(seed_cv)
  acc_CA   <- cv_grouped_K(tbv, list(KA_C),       cross_id, Kfold=Kfold, seed=seed_cv)
  acc_CAD  <- cv_grouped_K(tbv, list(KA_C, KD_C), cross_id, Kfold=Kfold, seed=seed_cv)
  
  acc_BA   <- cv_grouped_K(tbv, list(KA_B),       cross_id, Kfold=Kfold, seed=seed_cv)
  acc_BAD  <- cv_grouped_K(tbv, list(KA_B, KD_B), cross_id, Kfold=Kfold, seed=seed_cv)
  
  ## optional: variance components in tripPop (under SP4 trait)
  Va <- as.numeric(varA(tripPop, simParam = SP4))
  Vd <- as.numeric(varD(tripPop, simParam = SP4))
  
  data.frame(
    meanDD = meanDD,
    varDD  = varDD,
    h2_base = h2_base,
    Ve_fixed = Ve_fixed,
    n_train = n_cross * n_prog,   # total triploids
    n_cross = n_cross,
    n_prog  = n_prog,
    quadProb = quadProb,
    Va_trip = Va,
    Vd_trip = Vd,
    C_A   = mean(acc_CA,  na.rm=TRUE),
    C_AD  = mean(acc_CAD, na.rm=TRUE),
    Bmid_A  = mean(acc_BA,  na.rm=TRUE),
    Bmid_AD = mean(acc_BAD, na.rm=TRUE),
    gap_A   = mean(acc_CA,  na.rm=TRUE) - mean(acc_BA,  na.rm=TRUE),
    gap_AD  = mean(acc_CAD, na.rm=TRUE) - mean(acc_BAD, na.rm=TRUE)
  )
}

## ----------------------------
## 2) Run a small grid (edit these ranges)
## ----------------------------
## ----------------------------
## 2) Run a small grid (edit these ranges)
## ----------------------------
grid <- CJ(
  meanDD  = c(0.2, 0.6, 1.0),     # dominance degree
  varDD   = c(0.1, 0.5),          # dominance-degree variability
  h2_base = c(0.05, 0.12, 0.30),  # genetic signal
  n_cross = c(150, 400),          # training size lever #1
  n_prog  = c(2, 5)               # training size lever #2
)

## (Optional) keep quadProb fixed for now
grid[, quadProb := 0.3]

## ---------------------------------------------------------
## Run 50 replicates for every grid scenario
## ---------------------------------------------------------

n_reps <- 50

res <- rbindlist(
  lapply(seq_len(nrow(grid)), function(i){
    
    g <- grid[i]
    
    rbindlist(
      lapply(seq_len(n_reps), function(rep){
        
        cat(
          "Scenario", i, "of", nrow(grid),
          "| Replicate", rep, "of", n_reps, "\n"
        )
        
        out <- run_scenario(
          meanDD = g$meanDD,
          varDD  = g$varDD,
          h2_base = g$h2_base,
          n_cross = g$n_cross,
          n_prog  = g$n_prog,
          quadProb = g$quadProb,
          
          # make each replicate reproducible but different
          seed_main  = 29121983 + i*1000 + rep,
          seed_pheno = 111 + i*100 + rep,
          seed_cross = 9001 + i*10 + rep,
          seed_cv    = 123
        )
        
        # record replicate number
        out$rep <- rep
        
        out
      })
    )
  }),
  use.names = TRUE,
  fill = TRUE
)

print(res)

## save for plotting later
fwrite(res, "results_dosage_resolution_sensitivity_50reps.csv")

setwd("/Users/vmbl/Documents/Kajal_Simulation_2/results")
