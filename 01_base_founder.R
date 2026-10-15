### ============================
### CHUNK 0 — Setup + Base founder haplotypes + SP2 (diploid) with “locked” architecture
### ============================

suppressPackageStartupMessages({
  if (!requireNamespace("AlphaSimR", quietly = TRUE)) install.packages("AlphaSimR")
  library(AlphaSimR)
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
  library(data.table)
})

set.seed(29121983)

nChr        <- 11
segSites    <- 3500
nQTL        <- 495
chip.size.1 <- 3003


use_gamma    <- TRUE          
varA_target  <- 1             
meanDD_base  <- 0.6
varDD_base   <- 0.5         


h2_base   <- 0.12             
Ve_fixed  <- varA_target * (1 - h2_base) / h2_base  

cat(sprintf("Fixed Ve = %.4f (derived from h2_base=%.2f and varA_target=%.1f)\n",
            Ve_fixed, h2_base, varA_target))

qtl_per_chr <- nQTL / nChr        
snp_per_chr <- chip.size.1 / nChr 
stopifnot(qtl_per_chr == floor(qtl_per_chr),
          snp_per_chr == floor(snp_per_chr))

segSites_per_chr <- as.integer(segSites / nChr)  
if (segSites_per_chr * nChr != segSites) {
  message(sprintf("Note: using %d seg sites/chr (total %d), not %d.",
                  segSites_per_chr, segSites_per_chr*nChr, segSites))
}

history2 <- quickHaplo(
  nInd     = 1000,
  nChr     = nChr,
  segSites = segSites_per_chr
)


SP2 <- SimParam$new(history2)

SP2$restrSegSites(qtl_per_chr, snp_per_chr)

SP2$addTraitAD(
  nQtlPerChr = qtl_per_chr,
  mean       = 0,
  var        = varA_target,
  meanDD     = meanDD_base,
  varDD      = varDD_base,
  useVarA    = TRUE,
  gamma      = use_gamma
)

SP2$setVarE(varE = Ve_fixed)

# SNP chip
SP2$addSnpChip(snp_per_chr)

cat("CHUNK 0 OK\n")
cat("  nChr =", nChr,
    "| QTL/chr =", qtl_per_chr,
    "| SNP/chr =", snp_per_chr,
    "| segSites/chr =", segSites_per_chr, "\n")
cat("  SP2 fixed Ve =", round(Ve_fixed, 4),
    "| gamma =", use_gamma,
    "| meanDD =", meanDD_base, "\n")


### ============================
### CHUNK 1 — Create diploid founders
### ============================
firstPop <- newPop(rawPop = history2, isDH = FALSE, simParam = SP2)
firstPop@id <- paste0("A0_", sprintf("%03d", seq_len(nInd(firstPop))))


set.seed(111)
firstPop <- setPheno(firstPop, simParam = SP2)


cat("CHUNK 1 OK\n")
cat("  firstPop nInd =", nInd(firstPop), "| ploidy =", firstPop@ploidy, "\n")

Va <- as.numeric(varA(firstPop, simParam = SP2))
Vd <- as.numeric(varD(firstPop, simParam = SP2))
Vg <- Va + Vd               
Vp <- Vg + Ve_fixed           
h2_real <- Va / Vp

cat(sprintf("  VarA=%.4f | VarD=%.4f | VarG=%.4f | VarP=%.4f | h2(realized)=%.4f\n",
            Va, Vd, Vg, Vp, h2_real))


n.parents <- 96
set.seed(111)
firstPop <- setPheno(firstPop, simParam = SP2)

firstPop_sel <- selectInd(firstPop,
                          nInd = n.parents, trait = 1,
                          use = "pheno", gender = "B",
                          selectTop = TRUE,
                          returnPop = TRUE, simParam = SP2)
firstPop_sel@id <- paste0("A0_", sprintf("%03d", seq_len(nInd(firstPop_sel))))

cat("CHUNK 1.5 OK\n")
cat("  firstPop_sel nInd =", nInd(firstPop_sel), "| ploidy =", firstPop_sel@ploidy, "\n")

