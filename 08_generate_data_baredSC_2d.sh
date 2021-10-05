cd /scratch/ldelisle/LopezDelisle2021
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleEtAl2021/
mkdir -p output_baredSC_2d/mcmc/

conda activate baredSC

# First we generate data using the Nis from nih3t3

# 2d gauss:
python ${gitHubDirectory}/scripts/generate_2dgauss_data_columns.py --input output1/nih3t3_nRNA.txt --startingSeed 1 \
  --colnames "1_0.25_0.25_0.25_0.25_0.5" "1_0.25_0.25_0.25_0.25_-0.5" \
  "1_0.25_0.25_0.25_0.25_0" "1_0.25_0.25_0.15_0.25_0.5" \
  "1_0.25_0.25_0.15_0.25_-0.5" "1_0.25_0.25_0.15_0.25_0" \
  --output ${gitHubDirectory}/tables/nih3t3_generated_2d_1.txt

# Switching between 2 gaussians
python ${gitHubDirectory}/scripts/generate_paired_data_columns.py --input output1/nih3t3_nRNA.txt --startingSeed 1 \
  --colnamex1 "gauss_1_0.375_0.125" \
  --colnamex2 "gauss_1_1_0.1" \
  --colnamey1 "gauss_1_0.375_0.125" \
  --colnamey2 "gauss_1_1_0.1" \
  --props4groups "0.5_0_0_0.5" \
  "0_0.5_0.5_0" \
  --output ${gitHubDirectory}/tables/nih3t3_generated_2d_2.txt

# The only thing changing is the proportion of 0
python ${gitHubDirectory}/scripts/generate_paired_data_columns.py --input output1/nih3t3_nRNA.txt --startingSeed 2 \
  --colnamex1 "gauss_0_0.75_0.25" \
  --colnamex2 "gauss_1_0.75_0.25" \
  --colnamey1 "gauss_0_0.75_0.25" \
  --colnamey2 "gauss_1_0.75_0.25" \
  --props4groups "0.25_0.25_0.25_0.25" \
  "0.1_0.4_0.4_0.1" \
  "0_0.5_0.5_0" \
  --output ${gitHubDirectory}/tables/nih3t3_generated_2d_3.txt

conda deactivate

# Then we generate the table for parallel mcmc2d:
pathForTable="${gitHubDirectory}/tables/generated_table_2d.txt"
if [ -e $pathForTable ]; then
  rm $pathForTable
fi
nsim=0
table="${gitHubDirectory}/tables/nih3t3_generated_2d_1.txt"
my_simulations=$(head -n 1 ${table} | tr "\t" "\n" | grep -v expression | grep _x | sed "s/_x//g")
for simulation in $my_simulations; do
  echo -e "$table\t${simulation}_x\t${simulation}_y\t3\t3" >> $pathForTable
  echo -e "$table\t${simulation}_x\t${simulation}_y\t3\t3\tgroup" >> $pathForTable
  nsim=$((nsim + 2))
done

table="${gitHubDirectory}/tables/nih3t3_generated_2d_2.txt"
my_simulations_prop=$(head -n 1 ${table} | tr "\t" "\n" | grep -v expression | grep _x | sed "s/_x//g")
for simulation in $my_simulations_prop; do
  echo -e "$table\t${simulation}_x\t${simulation}_y\t3\t3" >> $pathForTable
  echo -e "$table\t${simulation}_x\t${simulation}_y\t3\t3\tgroup" >> $pathForTable
  nsim=$((nsim + 2))
done

table="${gitHubDirectory}/tables/nih3t3_generated_2d_3.txt"
my_simulations_prop=$(head -n 1 ${table} | tr "\t" "\n" | grep -v expression | grep _x | sed "s/_x//g")
for simulation in $my_simulations_prop; do
  echo -e "$table\t${simulation}_x\t${simulation}_y\t3\t3" >> $pathForTable
  echo -e "$table\t${simulation}_x\t${simulation}_y\t3\t3\tgroup" >> $pathForTable
  nsim=$((nsim + 2))
done

# I launch the parallel mcmc:
sbatch --array 1-${nsim} --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_2d.sh ${pathForTable} $PWD/output_baredSC_2d/mcmc/
