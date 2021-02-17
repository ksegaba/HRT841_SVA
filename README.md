# HRT841_SVA

File Descriptions:
- SVA_Entire_Dataset.Rmd: Surrogate Variable Analysis performed on the entire RNA-seq expression matrix (TPM units).
- CSS844_SVA_Group_PlottingCleanedMatrices.ipynb: Plotting expression matrices (TPM) that have been 'cleaned' of Surrogate Variable effects (generated from mdr_try_sva_v13.Rmd) using PCA, MDS, and tSNE.
- rPCA.R: Outlier removal using Robust Principle Component Analysis from the entire RNA-seq expression matrix before cleaning.
- fetch_bioproject.sh: Fetch corresponding BioProject number from an NCBI SRA Accession using Entrez Direct E-Utilities on the command line.

File links: 
Post-SVA and post-outlier-removal plots + datasets: https://drive.google.com/file/d/1qge-x1WB6csUD74AyMZmL9STc-BERZI4/view?usp=sharing
Outlier removal plots: https://drive.google.com/file/d/1LGS-RH6ApqrK1c8HI6mEbjKxHcRLrY2D/view?usp=sharing
