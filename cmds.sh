#Empty Drops
#sbatch Run_EmptyDrops.sh /cluster/projects/macparland/TA/Single_Nuc/Healthy/C58/ C58_TST
#sbatch Run_EmptyDrops.sh /cluster/projects/macparland/TA/Single_Nuc/Healthy/C70_NUC_TST_3pr_v3/raw_feature_bc_matrix/ C70_TST
sbatch Run_EmptyDrops.sh /cluster/projects/macparland/TA/Single_Nuc/Healthy/C72_NUC_TST_3pr_v3/raw_feature_bc_matrix/ C72_TST
#sbatch Run_EmptyDrops.sh /cluster/projects/macparland/TA/Single_Nuc/Healthy/SnRNA-seq_Data_C41/TST/ C41_TST
#sbatch Run_EmptyDrops.sh /cluster/projects/macparland/TA/Single_Nuc/Healthy/SnRNA-seq_Data_C41/NST/ C41_NST
#sbatch Run_EmptyDrops.sh /cluster/projects/macparland/TA/Single_Nuc/Healthy/SnRNA-seq_Data_C41/CST/ C41_CST
#sbatch Run_EmptyDrops.sh /cluster/projects/macparland/TA/Single_Nuc/PSC/PSC004S2 PSC_TST

#Pipeline

sbatch Run_Pipeline.sh /cluster/projects/macparland/TA/Single_Nuc/EmptyDrops/C41_CST_emptyDrops_table.rds
sbatch Run_Pipeline.sh /cluster/projects/macparland/TA/Single_Nuc/EmptyDrops/C41_NST_emptyDrops_table.rds
sbatch Run_Pipeline.sh /cluster/projects/macparland/TA/Single_Nuc/EmptyDrops/C41_TST_emptyDrops_table.rds
sbatch Run_Pipeline.sh /cluster/projects/macparland/TA/Single_Nuc/EmptyDrops/C58_TST_emptyDrops_table.rds
sbatch Run_Pipeline.sh /cluster/projects/macparland/TA/Single_Nuc/EmptyDrops/C70_TST_emptyDrops_table.rds
sbatch Run_Pipeline.sh /cluster/projects/macparland/TA/Single_Nuc/EmptyDrops/PSC_TST_emptyDrops_table.rds
