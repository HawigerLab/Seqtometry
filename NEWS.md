# Seqtometry 0.99.0

Update for potential Bioconductor submission

- Migrated quickstart from README.md to vignettes
- Updated function documentation with examples

# Seqtometry 0.2.0

Update in association with STAR Protocols paper

- Revised score function (improved performance by moving additional numerical code to Rcpp)
- Added impute function (R port of Python's magic-impute)

# Seqtometry 0.1.0

Initial release in association with Cell Systems paper

- Package provides function for signature (gene set) scoring tailored for single cell data
- Imputation relegated to magic-impute Python package (accessible from R via reticulate)
