### ============================
### CHUNK 2 — Build 4x founders + SP4 (LOCKED architecture from SP2)
### ============================

## --- requires from earlier chunks ---
stopifnot(exists("SP2"), exists("firstPop"))
stopifnot(exists("nChr"), exists("qtl_per_chr"), exists("snp_per_chr"), exists("Ve_fixed"))

genMap2 <- NULL
genMap2 <- tryCatch(SP2$genMap,  error = function(e) NULL)
if (is.null(genMap2)) genMap2 <- tryCatch(SP2$genoMap, error = function(e) NULL)
if (is.null(genMap2)) genMap2 <- tryCatch(history2@genMap, error = function(e) NULL)
stopifnot(!is.null(genMap2))

## Double genome to make 4x individuals (temporary)
firstPop4_tmp <- doubleGenome(firstPop_sel, simParam = SP2)  

## Convert to a proper MapPop 
H4 <- pullSegSiteHaplo(firstPop4_tmp, simParam = SP2)    
stopifnot(nrow(H4) == 4L * nInd(firstPop4_tmp))

cols_per_chr <- ncol(H4) / nChr
stopifnot(cols_per_chr == floor(cols_per_chr))

H4_list <- lapply(seq_len(nChr), function(i){
  a <- (i - 1) * cols_per_chr + 1
  b <- i * cols_per_chr
  H4[, a:b, drop = FALSE]
})

founders4 <- newMapPop(
  genMap     = genMap2,    
  haplotypes = H4_list,
  ploidy     = 4L
)


SP4 <- SimParam$new(founders4)

# Restrict eligible loci 
SP4$restrSegSites(qtl_per_chr, snp_per_chr)

## LOCK the exact same trait architecture from SP2 into SP4
SP4$manAddTrait(SP2$traits[[1]])

## Fix the same environmental variance (Ve)
SP4$setVarE(varE = Ve_fixed)

## Add SNP chip (positional arg)
SP4$addSnpChip(snp_per_chr)

## Make the tetraploid population under SP4
tetraPop <- newPop(rawPop = founders4, simParam = SP4)
tetraPop@id <- paste0("T0_", sprintf("%03d", seq_len(nInd(tetraPop))))

cat("CHUNK 2 OK\n")
cat("  tetraPop nInd =", nInd(tetraPop), "| ploidy =", tetraPop@ploidy, "\n")

## --- checks: trait locked? ---
cat("Same loci positions?  SP2 vs SP4:",
    isTRUE(all.equal(SP2$traits[[1]]@lociLoc, SP4$traits[[1]]@lociLoc)), "\n")
cat("Same additive effects? SP2 vs SP4:",
    isTRUE(all.equal(SP2$traits[[1]]@addEff,  SP4$traits[[1]]@addEff)), "\n")
cat("Same dominance effects? SP2 vs SP4:",
    isTRUE(all.equal(SP2$traits[[1]]@domEff,  SP4$traits[[1]]@domEff)), "\n")

## --- check Ve stored (some versions expose SP$varE) ---
cat("Fixed Ve check: SP2 varE =", round(SP2$varE, 4),
    "| SP4 varE =", round(SP4$varE, 4), "\n")

