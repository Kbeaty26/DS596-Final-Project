#!/bin/bash
#$ -P ds596
#$ -N cm_replication
#$ -l h_rt=01:00:00
#$ -pe omp 2
#$ -l mem_per_core=8G
#$ -o /projectnb/ds596/students/aliviap/DS596-Final-Project/extension/logs/replication.log
#$ -e /projectnb/ds596/students/aliviap/DS596-Final-Project/extension/logs/replication.err
#$ -j y

set -e
echo "Job started: $(date)"
echo "Host: $(hostname)"

module load R/4.5.2 2>/dev/null || export PATH="/share/pkg.8/r/4.5.2/install/bin:$PATH"
export R_LIBS_USER="/projectnb/ds596/students/aliviap/R_libs_4.5.2"

mkdir -p /projectnb/ds596/students/aliviap/DS596-Final-Project/extension/figures/replication
mkdir -p /projectnb/ds596/students/aliviap/DS596-Final-Project/extension/logs

Rscript /projectnb/ds596/students/aliviap/DS596-Final-Project/replication/replicate_figures.R

echo "Job finished: $(date)"