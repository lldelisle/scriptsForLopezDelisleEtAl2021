cd /scratch/ldelisle/LopezDelisle2021
module purge
module load gcc/7.4.0  openblas/0.3.6-openmp
module load r/3.6.0
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleEtAl2021/
Rscript ${gitHubDirectory}/scripts/plot_rouco_FACS.R ./ output_Pitx1/ output_Pitx1/Pitx1
#     mean1    sigma1 
# 1.2373491 0.1633461 
#        A1     mean1    sigma1        A2     mean2    sigma2     mean3    sigma3 
# 0.1486929 1.4158724 0.2514536 0.2321138 2.0547271 0.2005978 2.6341209 0.3820301 
#        A1     mean1    sigma1        A2     mean2    sigma2     mean3    sigma3 
# 0.2410368 1.3902565 0.1793438 0.1652642 1.7558717 0.1400436 2.2841845 0.4096456 


Rscript ${gitHubDirectory}/scripts/plot_1d_rouco_combinedModels.R ./ output_Pitx1/mcmc/ "${gitHubDirectory}/tables/rouco_table_1d.txt" output_Pitx1/Pitx1

Rscript ${gitHubDirectory}/scripts/plot_1d_rouco_combinedModels_log.R ./ output_Pitx1_log_15/mcmc/ "${gitHubDirectory}/tables/rouco_table_1d_log_15.txt" output_Pitx1_log_15/Pitx1
