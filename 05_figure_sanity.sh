cd /scratch/ldelisle/LopezDelisle2021
module purge
module load gcc/7.4.0  openblas/0.3.6-openmp
module load r/3.6.0
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleEtAl2021/
Rscript ${gitHubDirectory}/scripts/plot_1d_generated_1-4VS1gaussVSsanity.R ./ output_Sanity/mcmc/ "${gitHubDirectory}/tables/generated_table_1d_log.txt" output_Sanity/sanity output_Sanity/Sanity/
