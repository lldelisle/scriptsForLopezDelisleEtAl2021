cd /scratch/ldelisle/LopezDelisle2021
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleEtAl2021/
mkdir -p output_Sanity_2d/mcmc/

conda activate baredSC

# First we generate data using the Nis from nih3t3

# 2d gauss:
python ${gitHubDirectory}/scripts/generate_2dgauss_data_columns.py --input output1/nih3t3_nRNA.txt --startingSeed 1 \
  --colnames "1_-9_-9_0.5_0.5_0.5" "1_-9_-9_0.5_0.5_-0.5" \
  "1_-9_-9_0.5_0.5_0" \
  --xscale log --outputSanity ${gitHubDirectory}/tables/nih3t3_generated_2d_log_Sanity_1.txt \
  > ${gitHubDirectory}/tables/nih3t3_generated_2d_log_1.txt

# Switching between 2 gaussians
python ${gitHubDirectory}/scripts/generate_paired_data_columns.py --input output1/nih3t3_nRNA.txt --startingSeed 1 \
  --colnamex1 "gauss_1_-9_0.3" \
  --colnamex2 "gauss_1_-7.5_0.12" \
  --colnamey1 "gauss_1_-9_0.3" \
  --colnamey2 "gauss_1_-7.5_0.12" \
  --props4groups "0.5_0_0_0.5" \
  "0_0.5_0.5_0" \
  --xscale log --outputSanity ${gitHubDirectory}/tables/nih3t3_generated_2d_log_Sanity_2.txt \
  > ${gitHubDirectory}/tables/nih3t3_generated_2d_log_2.txt

# Switching between 2 gaussians
python ${gitHubDirectory}/scripts/generate_paired_data_columns.py --input output1/nih3t3_nRNA.txt --startingSeed 1 \
  --colnamex1 "gauss_1_-9_0.3" \
  --colnamex2 "gauss_1_-8_0.12" \
  --colnamey1 "gauss_1_-9_0.3" \
  --colnamey2 "gauss_1_-8_0.12" \
  --props4groups "0.5_0_0_0.5" \
  "0_0.5_0.5_0" \
  --xscale log --outputSanity ${gitHubDirectory}/tables/nih3t3_generated_2d_log_Sanity_3.txt \
  > ${gitHubDirectory}/tables/nih3t3_generated_2d_log_3.txt

conda deactivate

# Then we generate the table for parallel mcmc:
pathForTable="${gitHubDirectory}/tables/generated_table_2d_log.txt"
if [ -e $pathForTable ]; then
  rm $pathForTable
fi
nsim=0
xmin=-12
xmax=-6

table=${gitHubDirectory}/tables/nih3t3_generated_2d_log_1.txt
my_simulations=$(head -n 1 ${table} | tr "\t" "\n" | grep -v expression | grep _x | sed "s/_x//g")
for simulation in $my_simulations; do
  echo -e "${table}\t${simulation}_x\t${simulation}_y\t${xmin}\t${xmax}" >> $pathForTable
  nsim=$((nsim + 1))
done

table=${gitHubDirectory}/tables/nih3t3_generated_2d_log_2.txt
my_simulations_prop=$(head -n 1 ${table} | tr "\t" "\n" | grep -v expression | grep _x | sed "s/_x//g")
for simulation in $my_simulations_prop; do
  echo -e "${table}\t${simulation}_x\t${simulation}_y\t${xmin}\t${xmax}" >> $pathForTable
  nsim=$((nsim + 1))
done

table=${gitHubDirectory}/tables/nih3t3_generated_2d_log_3.txt
my_simulations_prop=$(head -n 1 ${table} | tr "\t" "\n" | grep -v expression | grep _x | sed "s/_x//g")
for simulation in $my_simulations_prop; do
  echo -e "${table}\t${simulation}_x\t${simulation}_y\t${xmin}\t${xmax}" >> $pathForTable
  nsim=$((nsim + 1))
done

# I launch the parallel mcmc2d in log scale:
sbatch --array 1-${nsim} --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_2d_log.sh ${pathForTable} $PWD/output_Sanity_2d/mcmc/

# I run sanity
mkdir -p output_Sanity_2d/Sanity_1/
Sanity/bin/Sanity -f ${gitHubDirectory}/tables/nih3t3_generated_2d_log_Sanity_1.txt -d output_Sanity_2d/Sanity_1/ -e true

mkdir -p output_Sanity_2d/Sanity_2/
Sanity/bin/Sanity -f ${gitHubDirectory}/tables/nih3t3_generated_2d_log_Sanity_2.txt -d output_Sanity_2d/Sanity_2/ -e true

mkdir -p output_Sanity_2d/Sanity_3/
Sanity/bin/Sanity -f ${gitHubDirectory}/tables/nih3t3_generated_2d_log_Sanity_3.txt -d output_Sanity_2d/Sanity_3/ -e true
