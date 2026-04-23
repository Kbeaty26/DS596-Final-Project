#!/bin/bash
#$ -P ds596
#$ -N cm_extension
#$ -l h_rt=02:00:00
#$ -pe omp 4
#$ -l mem_per_core=8G
#$ -o /projectnb/ds596/students/aliviap/extension/logs/extension.log
#$ -e /projectnb/ds596/students/aliviap/extension/logs/extension.err
#$ -j y

set -e

echo "Job started: $(date)"
echo "Host: $(hostname)"

# Load R
module load R/4.5.2 2>/dev/null || \
  export PATH="/share/pkg.8/r/4.5.2/install/bin:$PATH"

export R_LIBS_USER="/projectnb/ds596/students/aliviap/R_libs_4.5.2"

# Make sure log directory exists
mkdir -p /projectnb/ds596/students/aliviap/extension/logs
mkdir -p /projectnb/ds596/students/aliviap/extension/figures

echo "Running full extension demo (Baron -> Muraro)..."
Rscript /projectnb/ds596/students/aliviap/extension/run_extension_demo.R

echo "Job finished: $(date)"
