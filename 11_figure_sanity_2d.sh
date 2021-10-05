cd /scratch/ldelisle/LopezDelisle2021
module purge
module load gcc/7.4.0  openblas/0.3.6-openmp
module load r/3.6.0
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleEtAl2021/
Rscript ${gitHubDirectory}/scripts/plot_2d_generated_combinedModels_sanity.R ./ output_Sanity_2d/mcmc/ "${gitHubDirectory}/tables/generated_table_2d_log.txt" output_Sanity_2d/generated output_Sanity_2d/
