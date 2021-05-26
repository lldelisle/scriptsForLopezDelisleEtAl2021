cd /scratch/ldelisle/LopezDelisle2021
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleEtAl2021/
mkdir -p output_baredSC_1d/mcmc/

conda activate baredSC

# First we generate data using the Nis from nih3t3
python ${gitHubDirectory}/scripts/generate_data_columns.py --input output1/nih3t3_nRNA.txt --startingSeed 1 \
  --colnames "gauss_1_0.5_0.5" "gauss_1_0.75_0.25" "gauss_1_1_0.2" "gauss_1_1.5_0.5" \
  "uniform_1_0_1" "uniform_1_0_2" \
  "gauss_0.5_0.5_0.15_gauss_0.5_1.5_0.15" "gauss_0.5_0.75_0.25_gauss_0.5_2_0.2" \
  "gauss_0.5_1_0.5_gauss_0.5_2.5_0.15" \
  "gauss_0.3_0.5_0.15_gauss_0.4_1.5_0.15_gauss_0.3_2.5_0.15" \
  "gauss_0.3_0.5_0.2_gauss_0.3_1.25_0.2_gauss_0.4_2_0.2" \
  "gauss_0.25_0.75_0.25" "gauss_0.25_1.5_0.25" \
  "gauss_0.18_0.9_0.12" \
  > ${gitHubDirectory}/tables/nih3t3_generated_1d.txt

conda deactivate

# Then we generate the table for parallel mcmc:
pathForTable="${gitHubDirectory}/tables/generated_table_1d.txt"
my_simulations=$(head -n 1 ${gitHubDirectory}/tables/nih3t3_generated_1d.txt | tr "\t" "\n" | grep -v expression | grep -v group | grep -v nCount_RNA)
if [ -e $pathForTable ]; then
  rm $pathForTable
fi
nsim=0
for simulation in $my_simulations; do
  if [[ $simulation = *"2.5"* ]]; then
    xmax=4
  else
    xmax=3
  fi
  echo -e "${gitHubDirectory}/tables/nih3t3_generated_1d.txt\t${simulation}\t${xmax}" >> $pathForTable
  echo -e "${gitHubDirectory}/tables/nih3t3_generated_1d.txt\t${simulation}\t${xmax}\tgroup" >> $pathForTable
  nsim=$((nsim + 2))
done

# I launch the parallel mcmc:
sbatch --array 1-${nsim} --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d.sh ${pathForTable} $PWD/output_baredSC_1d/mcmc/
