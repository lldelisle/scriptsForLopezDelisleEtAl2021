cd /scratch/ldelisle/LopezDelisle2021
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleEtAl2021/
mkdir -p output_Pitx1/mcmc/
cp ${gitHubDirectory}/RoucoEtAl/* output_Pitx1/
# Then we generate the table for parallel mcmc:
pathForTable="${gitHubDirectory}/tables/rouco_table_1d.txt"
echo -e "output_Pitx1/RoucoEtAl.txt\tPitx1\t3\tlimbtype" > $pathForTable

# I launch the parallel mcmc:
sbatch --array 1 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d.sh ${pathForTable} $PWD/output_Pitx1/mcmc/

pathForTable="${gitHubDirectory}/tables/rouco_table_1d_log_15.txt"
echo -e "output_Pitx1/RoucoEtAl.txt\tPitx1\t-15\t-5\tlimbtype" > $pathForTable

mkdir -p output_Pitx1_log_15/mcmc/
cp ${gitHubDirectory}/RoucoEtAl/* output_Pitx1_log_15/
# I launch the parallel mcmc:
sbatch --array 1 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_log.sh ${pathForTable} $PWD/output_Pitx1_log_15/mcmc/
