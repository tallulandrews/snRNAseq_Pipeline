#!/bin/bash
#SBATCH -t 6-18:00:00
#SBATCH --mem=60000M 
#SBATCH -J Liver2_DE
#SBATCH -p himem 
#SBATCH -N 1 
#SBATCH -c 1 
#SBATCH -o %x-%j.out


echo $1 $2 $3 $4 $5 $6
#cd /cluster/projects/macparland/TA/Single_Nuc/SingleNuc_vs_SingleCell/Systematic/
#/cluster/tools/software/R/3.5.0/Rscript /cluster/home/tandrews/scripts/snRNAseq_Pipeline/Cell_type_markers.R $1 $2 $3 "Empty_Manual" "donor"

cd /cluster/projects/macparland/TA/Single_Nuc/SingleNuc_vs_SingleCell/
/cluster/tools/software/R/3.5.0/Rscript /cluster/home/tandrews/scripts/snRNAseq_Pipeline/Cell_type_markers.R "Sc_vs_sn_integrated_harmony_plus_analysis_Anno.rds" $2 $3 $4 $5 $6

# file
# clustername
# projname
# cluster_col_name (Use_clusters)
# confounder_col_name (donor)
# assay type

