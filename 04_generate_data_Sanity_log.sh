cd /scratch/ldelisle/LopezDelisle2021
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleEtAl2021/
mkdir -p output_Sanity/mcmc/

conda activate baredSC

# First we generate data using the Nis from nih3t3
python ${gitHubDirectory}/scripts/generate_data_columns.py --input output1/nih3t3_nRNA.txt --startingSeed 1 \
  --colnames "gauss_1_-9_0.5" "gauss_0.5_-9.5_0.5_gauss_0.5_-6_0.5" \
  "gauss_0.75_-9_0.5" \
  "gauss_0.5_-9_0.3_gauss_0.5_-7.5_0.1" \
  "gauss_0.3_-9.5_0.15_gauss_0.4_-8_0.15_gauss_0.3_-6.5_0.15" \
  --xscale log --outputSanity ${gitHubDirectory}/tables/nih3t3_generated_1d_log_Sanity.txt \
  > ${gitHubDirectory}/tables/nih3t3_generated_1d_log.txt

conda deactivate

# Then we generate the table for parallel mcmc:
pathForTable="${gitHubDirectory}/tables/generated_table_1d_log.txt"
my_simulations=$(head -n 1 ${gitHubDirectory}/tables/nih3t3_generated_1d_log.txt | tr "\t" "\n" | grep -v expression | grep -v group | grep -v nCount_RNA)
if [ -e $pathForTable ]; then
  rm $pathForTable
fi
nsim=0
xmin=-12
for simulation in $my_simulations; do
  xmax=-6
  if [[ $simulation = *"-6.5"* ]]; then
    xmax=-5
  fi
  if [[ $simulation = *"-6_0.5"* ]]; then
    xmax=-4
  fi
  echo -e "${gitHubDirectory}/tables/nih3t3_generated_1d_log.txt\t${simulation}\t${xmin}\t${xmax}" >> $pathForTable
  nsim=$((nsim + 1))
done
# To compare baredSC 1 gauss with Sanity we extend the x range:
echo -e "${gitHubDirectory}/tables/nih3t3_generated_1d_log.txt\tgauss_0.5_-9.5_0.5_gauss_0.5_-6_0.5\t-15\t-3" >> $pathForTable
nsim=$((nsim + 1))

# I launch the parallel mcmc in log scale:
sbatch --array 1-${nsim} --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_log.sh ${pathForTable} $PWD/output_Sanity/mcmc/
# I relaunch single gauss with pretty bins:
sbatch --array 1-${nsim} --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_1d_log_singleGauss.sh ${pathForTable} $PWD/output_Sanity/mcmc/

# I clone sanity
git clone https://github.com/jmbreda/Sanity.git
cd Sanity/src
make
cd ../..
# I run sanity
mkdir -p output_Sanity/Sanity/
Sanity/bin/Sanity -f ${gitHubDirectory}/tables/nih3t3_generated_1d_log_Sanity.txt -d output_Sanity/Sanity/ -e true

