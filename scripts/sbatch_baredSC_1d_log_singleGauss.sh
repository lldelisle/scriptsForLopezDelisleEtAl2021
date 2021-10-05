#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G # 4G * number of max nnorm * number of groups
#SBATCH --cpus-per-task 12 # number of max nnorm * number of groups
#SBATCH --time 24:00:00
#SBATCH --job-name baredSC_1d

pathForTable=$1
outputFolder=$2

module purge

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh

conda activate baredSC

# I don't know why but it is quicker when you set that:
# You need to evaluate how many
# MCMC will run at the same time to be sure to be below
# the number of cores available.
export OMP_NUM_THREADS=4

input=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
gene=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')
xmin=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $3}')
xmax=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $4}')
group=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $5}')

if [ -z $outputFolder ]; then
  outputFolder=./
fi

mkdir -p $outputFolder
nnorm=1

# First run the MCMCs
if [ -z $group ]; then
  title=${nnorm}gauss_${gene}_${SLURM_ARRAY_TASK_ID}
  baredSC_1d --input ${input} \
    --xscale log --xmin ${xmin} --xmax ${xmax} \
    --output ${outputFolder}/${title} --geneColName ${gene} \
    --nnorm ${nnorm} --minNeff 200 \
    --title ${title} --figure ${outputFolder}/${title}.png \
    --logevidence ${outputFolder}/${title}_logevid.txt --nis 10000 &> ${outputFolder}/${title}.log
else
  groupVals=$(less ${input} | awk -v g=$group 'NR==1{for(i=1;i<=NF;i++){if($i == g){col=i}}}NR>1{a[$col] +=1}END{for (v in a){print v}}')
  for val in $groupVals; do
    title=${nnorm}gauss_${gene}_${group}${val}_${SLURM_ARRAY_TASK_ID}
    baredSC_1d --input ${input} \
      --metadata1ColName ${group} --metadata1Values ${val} \
      --xscale log --xmin ${xmin} --xmax ${xmax} \
      --output ${outputFolder}/${title} --geneColName ${gene} \
      --nnorm ${nnorm} --minNeff 200 \
      --title ${title} --figure ${outputFolder}/${title}.png \
      --logevidence ${outputFolder}/${title}_logevid.txt --nis 10000 &> ${outputFolder}/${title}.log &
  done
fi
wait
# Make pretty
if [ -z $group ]; then
  title=${nnorm}gauss_${gene}_${SLURM_ARRAY_TASK_ID}
  baredSC_1d --input ${input} \
    --xscale log --xmin ${xmin} --xmax ${xmax} \
    --output ${outputFolder}/${title} --geneColName ${gene} \
    --nnorm ${nnorm} --minNeff 200 --prettyBins 200 \
    --title ${title} --figure ${outputFolder}/${title}_pretty.png &> ${outputFolder}/${title}_pretty.log
else
  groupVals=$(less ${input} | awk -v g=$group 'NR==1{for(i=1;i<=NF;i++){if($i == g){col=i}}}NR>1{a[$col] +=1}END{for (v in a){print v}}')
  for val in $groupVals; do
    title=${nnorm}gauss_${gene}_${group}${val}_${SLURM_ARRAY_TASK_ID}
    baredSC_1d --input ${input} \
      --metadata1ColName ${group} --metadata1Values ${val} \
      --xscale log --xmin ${xmin} --xmax ${xmax} \
      --output ${outputFolder}/${title} --geneColName ${gene} \
      --nnorm ${nnorm} --minNeff 200 --prettyBins 200 \
      --title ${title} --figure ${outputFolder}/${title}_pretty.png &> ${outputFolder}/${title}_pretty.log &
  done
fi