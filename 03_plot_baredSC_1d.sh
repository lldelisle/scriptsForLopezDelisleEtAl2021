cd /scratch/ldelisle/LopezDelisle2021
module purge
module load gcc/7.4.0  openblas/0.3.6-openmp
module load r/3.6.0
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleEtAl2021/
Rscript ${gitHubDirectory}/scripts/plot_1d_generated_combinedModels.R ./ output_baredSC_1d/mcmc/ "${gitHubDirectory}/tables/generated_table_1d.txt" output_baredSC_1d/generated
