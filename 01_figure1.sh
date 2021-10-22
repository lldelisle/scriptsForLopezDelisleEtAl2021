cd /scratch/ldelisle/LopezDelisle2021
module purge
module load gcc/7.4.0  openblas/0.3.6-openmp
module load hdf5/1.10.5
module load r/3.6.0
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleEtAl2021/
mkdir -p output1
# This requires quite some memory
# For some configuration issues on the cluster, the first part until write.table(big.df) was performed with R 3.6.0
# While the figures were done on R 4.1.1 with ggplot2 3.3.5
Rscript ${gitHubDirectory}/scripts/check_Poisson_controlRNA.R $PWD/ output1
