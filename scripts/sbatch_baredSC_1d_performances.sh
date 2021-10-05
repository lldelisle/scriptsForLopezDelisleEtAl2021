#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --time 30:00:00
#SBATCH --job-name baredSC_1d_perfo

pathForTable=$1
outputFolder=$2
mingauss=$3
maxgauss=$4
seed=$5
if [ -z $mingauss ]; then
  mingauss=1
fi
if [ -z $maxgauss ]; then
  maxgauss=4
fi
if [ -z $seed ]; then
  seed=1
fi

module purge

source $(dirname $(dirname $(which conda)))/etc/profile.d/conda.sh

conda activate baredSC

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

input=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
gene=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')
xmax=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $3}')
group=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $4}')

if [ -z $outputFolder ]; then
  outputFolder=./
fi

mkdir -p $outputFolder

# To increase the writting speed a temp folder is created and then results are copied back
final_output=$outputFolder
outputFolder=$(mktemp -d)

# First run the MCMCs
if [ -z $group ]; then
  for (( nnorm=${mingauss}; nnorm<=${maxgauss}; nnorm++ )); do
    title=${nnorm}gauss_${gene}_${SLURM_ARRAY_TASK_ID}
    baredSC_1d --input ${input} \
      --xmax ${xmax} \
      --output ${outputFolder}/${title} --geneColName ${gene} \
      --nnorm ${nnorm} --minNeff 200 --seed ${seed} \
      --title ${title} --figure ${outputFolder}/${title}.png \
      --logevidence ${outputFolder}/${title}_logevid.txt --nis 10000 &> ${outputFolder}/${title}.log &
  done
else
  groupVals=$(less ${input} | awk -v g=$group 'NR==1{for(i=1;i<=NF;i++){if($i == g){col=i}}}NR>1{a[$col] +=1}END{for (v in a){print v}}')
  for val in $groupVals; do
    for (( nnorm=${mingauss}; nnorm<=${maxgauss}; nnorm++ )); do
      title=${nnorm}gauss_${gene}_${group}${val}_${SLURM_ARRAY_TASK_ID}
      baredSC_1d --input ${input} \
        --metadata1ColName ${group} --metadata1Values ${val} \
        --xmax ${xmax} \
        --output ${outputFolder}/${title} --geneColName ${gene} \
        --nnorm ${nnorm} --minNeff 200 --seed ${seed} \
        --title ${title} --figure ${outputFolder}/${title}.png \
        --logevidence ${outputFolder}/${title}_logevid.txt --nis 10000 &> ${outputFolder}/${title}.log &
    done
  done
fi
wait
if [ ! "$mingauss" = "$maxgauss" ]; then
  # Combine them:
  if [ -z $group ]; then
    outputs=""
    for (( nnorm=${mingauss}; nnorm<=${maxgauss}; nnorm++ )); do
      title=${nnorm}gauss_${gene}_${SLURM_ARRAY_TASK_ID}
      outputs="$outputs ${outputFolder}/${title}"
    done
    title=${mingauss}-${maxgauss}gauss_${gene}_${SLURM_ARRAY_TASK_ID}
    combineMultipleModels_1d --input ${input} \
      --xmax ${xmax} \
      --outputs ${outputs} --geneColName ${gene} \
      --title ${title} --figure ${outputFolder}/${title}.png \
      --nis 10000 &> ${outputFolder}/${title}.log &
  else
    groupVals=$(less ${input} | awk -v g=$group 'NR==1{for(i=1;i<=NF;i++){if($i == g){col=i}}}NR>1{a[$col] +=1}END{for (v in a){print v}}')
    for val in $groupVals; do
      outputs=""
      for (( nnorm=${mingauss}; nnorm<=${maxgauss}; nnorm++ )); do
        title=${nnorm}gauss_${gene}_${group}${val}_${SLURM_ARRAY_TASK_ID}
        outputs="$outputs ${outputFolder}/${title}"
      done
      title=${mingauss}-${maxgauss}gauss_${gene}_${group}${val}_${SLURM_ARRAY_TASK_ID}
      combineMultipleModels_1d --input ${input} \
        --metadata1ColName ${group} --metadata1Values ${val} \
        --xmax ${xmax} \
        --outputs ${outputs} --geneColName ${gene} \
        --title ${title} --figure ${outputFolder}/${title}.png \
        --nis 10000 &> ${outputFolder}/${title}.log &
    done
  fi
  wait
fi
# Copy back the data
rsync -r $outputFolder/ $final_output
