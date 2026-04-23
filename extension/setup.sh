#!/bin/bash
set -e

# Detect username and set personal R library path
USERNAME=$(whoami)
R_LIB="/projectnb/ds596/students/${USERNAME}/R_libs_4.5.2"

echo "=== CellMentor Extension Setup ==="
echo "Username:  ${USERNAME}"
echo "R library: ${R_LIB}"
echo ""

# Load R
module load R/4.5.2 2>/dev/null || {
  echo "Could not load R/4.5.2 via module. Trying path fallback..."
  export PATH="/share/pkg.8/r/4.5.2/install/bin:$PATH"
}

# Create personal library directory
mkdir -p "${R_LIB}"
export R_LIBS_USER="${R_LIB}"

echo "Installing dependencies (skipping already-installed packages)..."
echo ""

Rscript - <<REOF
.libPaths(c("${R_LIB}", .libPaths()))

cran_pkgs <- c(
  "Matrix", "ggplot2", "dplyr", "patchwork", "tidyr",
  "aricode", "pROC", "devtools", "BiocManager",
  "RMTstat", "skmeans", "MLmetrics", "lsa"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing: ", p)
    install.packages(p,
      lib   = "${R_LIB}",
      repos = "https://cloud.r-project.org",
      quiet = FALSE
    )
  } else {
    message("Already installed: ", p)
  }
}

# Seurat via CRAN
if (!requireNamespace("Seurat", quietly = TRUE)) {
  message("Installing Seurat...")
  install.packages("Seurat",
    lib   = "${R_LIB}",
    repos = "https://cloud.r-project.org"
  )
} else {
  message("Already installed: Seurat")
}

# CellMentor via Bioconductor (preferred) or GitHub fallback
if (!requireNamespace("CellMentor", quietly = TRUE)) {
  message("Installing CellMentor via Bioconductor...")
  ok <- tryCatch({
    BiocManager::install("CellMentor",
      lib    = "${R_LIB}",
      ask    = FALSE,
      update = FALSE
    )
    TRUE
  }, error = function(e) {
    message("Bioconductor install failed, trying GitHub...")
    FALSE
  })
  if (!ok) {
    devtools::install_github("petrenkokate/CellMentor",
      lib          = "${R_LIB}",
      dependencies = TRUE
    )
  }
} else {
  message("Already installed: CellMentor")
}

# Verify everything
cat("\n=== Verification ===\n")
pkgs <- c("CellMentor", "Seurat", "Matrix", "ggplot2", "dplyr",
          "patchwork", "aricode", "pROC", "RMTstat", "skmeans",
          "MLmetrics", "lsa")
all_ok <- TRUE
for (p in pkgs) {
  ok <- requireNamespace(p, quietly = TRUE)
  if (!ok) all_ok <- FALSE
  cat(sprintf("  %-20s %s\n", p, if (ok) "OK" else "MISSING"))
}
if (all_ok) {
  cat("\nAll packages installed successfully.\n")
} else {
  cat("\nSome packages are missing. Re-run setup.sh or install manually.\n")
  quit(status = 1)
}
REOF

echo ""
echo "=== Setup complete ==="
echo ""
echo "Your R library is at: ${R_LIB}"
echo ""
echo "Before running any scripts, either:"
echo "  export R_LIBS_USER=\"${R_LIB}\""
echo "or ensure the scripts reference this path (already done if you"
echo "used the automated sed commands in the README)."
