# gch1-variant-landscape
Computational characterisation of GCH1 missense variants
# GCH1 Variant Landscape

A computational analysis and comparison of GCH1 missense variants sourced from ClinVar.
This project aims to compare physicochemical and conservation-based features using parameters like BLOSUM62, Grantham distance, Δcharge and Δhydrophobicity in pathogenic 
and benign variants.  The pathogenic R198Q variant will then be compared in context of the broader variant landscape.

## Tools used
Python, pandas, scipy, scikit-learn, biopython, matplotlib, seaborn

## Data Sources
- ClinVar (https://www.ncbi.nlm.nih.gov/clinvar/)
- UniProt (https://www.uniprot.org/)

## Structure
- data/raw — raw downloaded files
- data/processed — cleaned and feature-annotated datasets
- notebooks — analysis notebooks
- scripts — reusable Python functions
- figures — output plots
- docs — written summaries and interpretation