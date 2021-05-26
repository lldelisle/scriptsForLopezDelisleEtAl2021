cd /scratch/ldelisle/LopezDelisle2021
gitHubDirectory=/home/ldelisle/softwares/scriptsForLopezDelisleEtAl2021/
mkdir -p output_Bolt/mcmc/
# Process data from Bolt et al.
wget https://github.com/lldelisle/scriptsForBoltEtAl2021/raw/45887c3e3c497a3c94e5f853c128fdbff60bc38f/scRNAseq/meta_data_scRNAseq.txt.gz
gunzip meta_data_scRNAseq.txt.gz
mv meta_data_scRNAseq.txt output_Bolt/
# Only keep distal clusters: 3 4 5 9 and wt genotype
cat output_Bolt/meta_data_scRNAseq.txt | awk -v cl="3,4,5,9" '
BEGIN{
  split(cl, a, ",")
  for(j in a){
    valid_cl[a[j]]=1
  }
}
NR==1{
  for(k=1;k<=NF;k++){
    if($k == "seurat_clusters"){
      # The table was exported with raw names:
      i=k+1
    }
    if($k == "genotype"){
      # The table was exported with raw names:
      j=k+1
    }
  }
  print
}
{
  if (($i in valid_cl) && ($j == "wt")){
      print
  }
}' > output_Bolt/distal_clusters.txt
# Then we generate the table for parallel mcmc:
pathForTable="${gitHubDirectory}/tables/bolt_table_2d.txt"
echo -e "output_Bolt/distal_clusters.txt\tHoxd13\tHoxa11\t3\t3\tseurat_clusters" > $pathForTable
# I launch the parallel mcmc:
sbatch --array 1 --chdir $PWD/ ${gitHubDirectory}/scripts/sbatch_baredSC_2d.sh ${pathForTable} $PWD/output_Bolt/mcmc/
