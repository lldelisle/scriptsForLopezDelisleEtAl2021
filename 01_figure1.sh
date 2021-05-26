cd /scratch/ldelisle/LopezDelisle2021
module purge
module load gcc/7.4.0  openblas/0.3.6-openmp
module load hdf5/1.10.5
module load r/3.6.0
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleEtAl2021/
mkdir -p output1
# This requires quite some memory
Rscript ${gitHubDirectory}/scripts/check_Poisson_controlRNA.R $PWD/ output1
