# =============================================================================
# run_pipeline.R
# Master pipeline that runs the full study end to end
#
# Run from Terminal (macOS):
#   caffeinate -i Rscript --vanilla code/run_pipeline.R 2>&1 | tee logs/master.log
#
# Expected wall time: 8-12 hours on a modern desktop
# =============================================================================

# OpenMP pinning to prevent macOS Tahoe deadlocks in AlphaSimR
Sys.setenv(OMP_NUM_THREADS = "2")
Sys.setenv(OMP_DYNAMIC = "FALSE")
Sys.setenv(OMP_NESTED = "FALSE")

# Timestamped logging helper
say <- function(msg) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))
  flush.console()
}

t_master <- Sys.time()
say("=== PIPELINE START ===")
say(sprintf("R version: %s", R.version$version.string))
say(sprintf("AlphaSimR: %s", packageVersion("AlphaSimR")))

# Create output directories
dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

# ---- Stage 1: diploid founders and SP2 ----
say(">>> Stage 1: 01_base_founder.R")
t0 <- Sys.time()
source("code/01_base_founder.R", echo = FALSE)
say(sprintf("Stage 1 done in %.1f min",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- Stage 2: tetraploid founders and SP4 ----
say(">>> Stage 2: 02_build_4x_founders.R")
t0 <- Sys.time()
source("code/02_build_4x_founders.R", echo = FALSE)
say(sprintf("Stage 2 done in %.1f min",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Save setup objects as a checkpoint
say("Saving setup objects to results/setup.rds")
saveRDS(list(
  SP2 = SP2, SP4 = SP4,
  firstPop = firstPop, firstPop_sel = firstPop_sel,
  tetraPop = tetraPop, Ve_fixed = Ve_fixed,
  nChr = nChr, qtl_per_chr = qtl_per_chr, snp_per_chr = snp_per_chr
), "results/setup.rds")

# ---- Stage 3: 50-rep sensitivity analysis (long) ----
say(">>> Stage 3: 03_sensitivity_50reps.R (~8 hours)")
t0 <- Sys.time()
source("code/03_sensitivity_50reps.R", echo = FALSE)
say(sprintf("Stage 3 done in %.1f min",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- Stage 4: 100-rep recurrent selection (long) ----
say(">>> Stage 4: run_RS_v2_100.R (~3 hours)")
t0 <- Sys.time()
source("code/run_RS_v2_100.R", echo = FALSE)
say(sprintf("Stage 4 done in %.1f min",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- Stage 5: statistical analysis of RS results ----
say(">>> Stage 5: analyze_RS.R")
t0 <- Sys.time()
source("code/analyze_RS.R", echo = FALSE)
say(sprintf("Stage 5 done in %.1f min",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- Stage 6: sensitivity figures ----
say(">>> Stage 6: Figures_50rep.R")
t0 <- Sys.time()
source("code/Figures_50rep.R", echo = FALSE)
say(sprintf("Stage 6 done in %.1f min",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

say(sprintf("=== PIPELINE DONE in %.1f min total ===",
            as.numeric(difftime(Sys.time(), t_master, units = "mins"))))