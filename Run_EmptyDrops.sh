#!/bin/bash
#SBATCH -t 6-18:00:00
#SBATCH --mem=80000M 
#SBATCH -J snEmpty
#SBATCH -p himem 
#SBATCH -N 1 
#SBATCH -c 1 
#SBATCH -o %x-%j.out



FOLDER=$1
PROJ=$2

cd /cluster/projects/macparland/TA/Single_Nuc

SCRIPT=/cluster/home/tandrews/scripts/snRNAseq_Pipeline/Recall_nuclei_with_EmptyDrops.R


echo $FOLDER
echo $PROJ

/cluster/tools/software/R/3.5.0/Rscript $SCRIPT $FOLDER $PROJ
