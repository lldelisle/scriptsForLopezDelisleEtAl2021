cd /scratch/ldelisle/LopezDelisle2021
module purge
module load gcc/7.4.0  openblas/0.3.6-openmp
module load r/3.6.0
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleEtAl2021/
Rscript ${gitHubDirectory}/scripts/plot_2d_bolt_combinedModels.R ./ output_Bolt/mcmc/ "${gitHubDirectory}/tables/bolt_table_2d.txt" output_Bolt/Bolt
