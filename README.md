# Seqtometry
This repository contains code for the Seqtometry scoring algorithm from the following publication:

Kousnetsov et al., Single-cell sequencing analysis within biologically relevant dimensions, Cell Systems (2023), https://doi.org/10.1016/j.cels.2023.12.005

For the accompanying GUI prototype, please see the [SeqtometryGUI repository](https://github.com/HawigerLab/SeqtometryGUI).

## Installation
This package can be installed using devtools.
```r
devtools::install_github("HawigerLab/Seqtometry")
```

## Usage notes
Note that the Seqtometry::score function assumes that the user has already performed quality control, normalization, and imputation
as described in the Methods section of the Seqtometry publication (linked above).
