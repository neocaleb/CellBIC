# CellBIC

CellBIC is a tool to cluster single cell transcriptomic data into top-down hierarchical clusters using bimodality in the gene expression distribution. CellBIC is implemented in MATLAB.

## Usage

1. run CellBIC_step1 with log-transformed count data and parameters. This function will return a top-down hierarchical clustering result.
2. run CellBIC_step2 with the returned variables from CellBIC_step1. THis function will return a clustering result with a given number of clusters.

## Example

Four example scripts are available with the corresponding single cell RNA sequencing data as follows:
1. RunCellBIC_Enge.m [1]
2. RunCellBIC_Treutlein.m [2]
3. RunCellBIC_Wang.m [3]
4. RunCellBIC_Zeisel.m [4]

## Citation

Manuscript for CellBIC is available from Nucleic Acids Research [5].

## References
1. Enge, M. et al. Single-Cell Analysis of Human Pancreas Reveals Transcriptional Signatures of Aging and Somatic Mutation Patterns. Cell 171, 321–330.e14 (2017).
2. Treutlein, B. et al. Reconstructing lineage hierarchies of the distal lung epithelium using single-cell RNA-seq. Nature 509, 371–5 (2014).
3. Wang, Y. J. et al. Single cell transcriptomics of the human endocrine pancreas. Diabetes (2016).
4. Zeisel,  a. et al. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. Science 347, 1138–42 (2015).
5. Kim, J. et al. CellBIC: bimodality-based top-down clustering of single-cell RNA sequencing data reveals hierarchical structure of the cell type. Nucelic Acids Research 46 (21), e124-e124 (2018).

