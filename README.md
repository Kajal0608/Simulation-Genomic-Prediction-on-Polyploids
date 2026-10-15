# Simulation study of allele dosage encoding in triploid genomic prediction

Companion code and data for:

Gupta K, Grover K, Pereira H, Fritsche-Neto R. 2026. The importance of allele dosage: a simulation study on genomic prediction and response to selection in triploid hybrid breeding. Theoretical and Applied Genetics (under review).

## Overview

This repository contains R scripts and output data used in the paper. The simulation compares two triploid genotype encodings for genomic prediction:

- Encoding C: proper four-state dosage (0, 1, 2, 3)
- Encoding B: diploidized three-state encoding (0, 1, 1, 2)

Two experiments are reproduced:

1. A 72-scenario factorial sensitivity analysis of prediction accuracy, replicated 50 times per scenario.
2. A 10-cycle recurrent selection program, replicated 100 times with paired seeds so that both encodings receive identical stochastic realizations.

## Requirements

- R 4.5.1 or later
- R packages: AlphaSimR, rrBLUP, data.table, ggplot2, patchwork, scales

Install with:

install.packages(c("AlphaSimR", "rrBLUP", "data.table", "ggplot2", "patchwork", "scales"))

## Reproducing the paper

To reproduce the full study:

Rscript run_pipeline.R

Expected wall time: 8 to 12 hours on a modern desktop.

## Contact

Kajal Gupta, corresponding author
Department of Horticultural Science, North Carolina State University

## License

MIT License.
