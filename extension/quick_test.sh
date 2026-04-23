#!/bin/bash
# quick_test.sh

export R_LIBS_USER="/projectnb/ds596/students/aliviap/R_libs_4.5.2"
module load R/4.5.2 2>/dev/null || \
  export PATH="/share/pkg.8/r/4.5.2/install/bin:$PATH"

mkdir -p /projectnb/ds596/students/aliviap/extension/figures
mkdir -p /projectnb/ds596/students/aliviap/extension/logs

Rscript /projectnb/ds596/students/aliviap/extension/run_extension_demo.R --subset \
  2>&1 | tee /projectnb/ds596/students/aliviap/extension/logs/quick_test.log

echo ""
echo "Check outputs:"
ls -lh /projectnb/ds596/students/aliviap/extension/figures/
