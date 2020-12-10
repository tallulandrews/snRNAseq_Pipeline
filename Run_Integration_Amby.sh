#!/bin/bash
#SBATCH -t 6-18:00:00
#SBATCH --mem=60000M 
#SBATCH -J integration_sn
#SBATCH -p himem 
#SBATCH -N 1 
#SBATCH -c 1 
#SBATCH -o %x-%j.out


cd /cluster/projects/macparland/TA/Single_Nuc/Amby
/cluster/tools/software/R/3.5.0/Rscript /cluster/home/tandrews/scripts/snRNAseq_Pipeline/Harmony_Integration_Amby.R
