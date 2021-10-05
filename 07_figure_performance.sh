cd /scratch/ldelisle/LopezDelisle2021
module purge
module load gcc/7.4.0  openblas/0.3.6-openmp
module load r/3.6.0
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleEtAl2021/

while read line; do
  jobid=$(echo $line | awk '{print $1}')
  outputfile=$(echo $line | awk '{print $2}')
  sacct -j ${jobid} --format=JobID%15,Elapsed,MaxVMSize,MaxRSS,State --noconvert --parsable2 > ${outputdir}/${outputfile}
done < ${gitHubDirectory}/tables/jobid_correspondance.txt

Rscript ${gitHubDirectory}/scripts/plot_benchmark.R ./ output_baredSC_1d_perf/ "${gitHubDirectory}/tables/generated_table_1d_perf.txt" output_baredSC_1d_perf/fig_benchmark

