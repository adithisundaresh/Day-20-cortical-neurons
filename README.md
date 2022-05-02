# Day-20-cortical-neurons
Master's thesis work
Processing of scRNA-seq data from day 20 _in vitro_ differentiated cortical neural progenitors

## Input
Read10X takes sample_feature_bc_matrix files from per_sample_outs directory of cell ranger multi output

## Files
1. [**Data_preprocessing.R**](Data_preprocessing.R) - merging, filtering and pre-processing steps
2. [**Annotation.R**](Annotation.R) - annotation with custom markers and fetal reference
3. [**Reproducibility_tests.R**](Reproducibility_tests.R) - evaluation of reproducibility between batches, cell lines and clones
