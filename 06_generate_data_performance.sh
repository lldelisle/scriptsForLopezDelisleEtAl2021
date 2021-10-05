cd /scratch/ldelisle/LopezDelisle2021
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleEtAl2021/
mkdir -p output_baredSC_1d_perf
outputdir=${PWD}/output_baredSC_1d_perf/

conda activate baredSC
mytmp=$(mktemp -d)
ns="1000 5000 25000 125000 625000 3125000"
for n in $ns; do
  python ${gitHubDirectory}/scripts/generate_nCounts.py --ncells $n > ${mytmp}/counts${n}.txt
  python ${gitHubDirectory}/scripts/generate_data_columns.py --input ${mytmp}/counts${n}.txt \
  --colnames "gauss_1_0.75_0.25" "gauss_0.5_0.5_0.15_gauss_0.5_1.5_0.15" \
  "gauss_0.3_0.5_0.15_gauss_0.4_1.5_0.15_gauss_0.3_2.5_0.15" \
  "gauss_0.25_0.75_0.25" "gauss_0.25_1.5_0.25" \
  > ${outputdir}/generated_1d_${n}.txt
done

conda deactivate

# Then we generate the table for parallel mcmc:
pathForTable="${gitHubDirectory}/tables/generated_table_1d_perf.txt"
my_simulations=$(head -n 1 ${outputdir}/generated_1d_1000.txt | tr "\t" "\n" | grep -v expression | grep -v group | grep -v nCount_RNA)
if [ -e $pathForTable ]; then
  rm $pathForTable
fi
nsim=0
for n in $ns; do
  for simulation in $my_simulations; do
    if [[ $simulation = *"2.5"* ]]; then
      xmax=4
    else
      xmax=3
    fi
    echo -e "${outputdir}/generated_1d_${n}.txt\t${simulation}\t${xmax}" >> $pathForTable
    nsim=$((nsim + 1))
  done
done
# I launch the parallel mcmc:
sbatch --array 1-${nsim} --ntasks 4 --cpus-per-task 4 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_default_4cores/
# job id: 1825496
sbatch --array 1-${nsim} --ntasks 4 --cpus-per-task 4 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_default_4cores_2/ 1 4 2
# job id: 1825703
sbatch --array 1-${nsim} --ntasks 4 --cpus-per-task 4 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_default_4cores_3/ 1 4 3
# job id: 1826276
#####
# Then only the expected model 4 cores:
sbatch --array 1,6,11,16,21,26 --ntasks 1 --cpus-per-task 4 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_4cores/ 1 1
# 1825528
sbatch --array 2,4,5,7,9,10,12,14,15,17,19,20,22,24,25,27,29,30 --ntasks 1 --cpus-per-task 4 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_4cores/ 2 2
# 1825529
sbatch --array 3,8,13,18,23,28 --ntasks 1 --cpus-per-task 4 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances48h.sh ${pathForTable} ${outputdir}/mcmc_expectedM_4cores/ 3 3
# 1825530
sbatch --array 1,6,11,16,21,26 --ntasks 1 --cpus-per-task 4 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_4cores/ 1 1 2
# 1825765
sbatch --array 2,4,5,7,9,10,12,14,15,17,19,20,22,24,25,27,29,30 --ntasks 1 --cpus-per-task 4 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_4cores/ 2 2 2
# 1825766
sbatch --array 3,8,13,18,23,28 --ntasks 1 --cpus-per-task 4 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_4cores/ 3 3 2
# 1825767
sbatch --array 1,6,11,16,21,26 --ntasks 1 --cpus-per-task 4 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_4cores_3/ 1 1 3
# 1825848
sbatch --array 2,4,5,7,9,10,12,14,15,17,19,20,22,24,25,27,29,30 --ntasks 1 --cpus-per-task 4 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_4cores_3/ 2 2 3
# 1825849
sbatch --array 3,8,13,18,23,28 --ntasks 1 --cpus-per-task 4 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_4cores_3/ 3 3 3
# 1825850
# Only the expected 1 core:
sbatch --array 1,6,11,16,21 --ntasks 1 --cpus-per-task 1 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_1core/ 1 1
# 1825907
sbatch --array 2,4,5,7,9,10,12,14,15,17,19,20,22,24,25 --ntasks 1 --cpus-per-task 1 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_1core/ 2 2
# 1825908
sbatch --array 3,8,13,18,23 --ntasks 1 --cpus-per-task 1 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_1core/ 3 3
# 1825909
sbatch --array 1,6,11,16,21 --ntasks 1 --cpus-per-task 1 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_1core_2/ 1 1 2
# 1825925
sbatch --array 2,4,5,7,9,10,12,14,15,17,19,20,22,24,25 --ntasks 1 --cpus-per-task 1 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_1core_2/ 2 2 2
# 1825926
sbatch --array 3,8,13,18,23 --ntasks 1 --cpus-per-task 1 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_1core_2/ 3 3 2
# 1825927
sbatch --array 1,6,11,16,21 --ntasks 1 --cpus-per-task 1 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_1core_3/ 1 1 3
# 1826175
sbatch --array 2,4,5,7,9,10,12,14,15,17,19,20,22,24,25 --ntasks 1 --cpus-per-task 1 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_1core_3/ 2 2 3
# 1826176
sbatch --array 3,8,13,18,23 --ntasks 1 --cpus-per-task 1 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_performances.sh ${pathForTable} ${outputdir}/mcmc_expectedM_1core_3/ 3 3 3
# 1826177
