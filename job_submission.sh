#!/bin/bash --login
#SBATCH --time=3:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=36G
#SBATCH --job-name run_sva
#SBATCH --output=%x_%j

cd /mnt/home/seguraab/HRT841/SVA_Analysis/ # path to working directory

module purge
module load GCC/6.4.0-2.28  OpenMPI/2.1.2  R/3.5.1-X11-20180131

# Get a list of SRAs with expression data
Rscript get_SRAs_with_expression.r Orthogroup_RNAseq_8-23-21.csv metadata_labels.csv

# Obtain Bioproject IDs for the SRAs
# edirect must be installed and in your PATH variable
# Installation: https://www.ncbi.nlm.nih.gov/books/NBK179288/
bash fetch_bioproject.sh
# Delete the extra bioproject IDs in the second column, if any, before the next step

# Run SVA
Rscript SVA_Entire_Dataset.r Orthogroup_RNAseq_8-23-21.csv metadata_labels.csv bioproject_ids.txt NEW_stress NEW_tissue NEW_family Orthogroup_RNAseq_8-23-21

scontrol show job $SLURM_JOB_ID # job information
