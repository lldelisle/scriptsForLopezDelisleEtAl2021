# On helvetios
mkdir /scratch/ldelisle/LopezDelisle2021
cd /scratch/ldelisle/LopezDelisle2021
# Download the h5ad from the figshare https://figshare.com/projects/Zero_inflation_in_negative_control_data/61292 kindly provided by Valentine Svensson
wget https://ndownloader.figshare.com/files/14634407 -O zheng_gemcode_control.h5ad
wget https://ndownloader.figshare.com/files/14634410 -O svensson_chromium_control.h5ad
wget https://ndownloader.figshare.com/files/14634416 -O klein_indrops_control_GSM1599501.h5ad
wget https://ndownloader.figshare.com/files/14634488 -O macosko_dropseq_control.h5ad
# Download data used in Svensson, V. Droplet scRNA-seq is not zero-inflated. Nat Biotechnol 38, 147–150 (2020). https://doi.org/10.1038/s41587-019-0379-5
wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/hgmm_5k_v3/hgmm_5k_v3_filtered_feature_bc_matrix.tar.gz
tar zxvmf hgmm_5k_v3_filtered_feature_bc_matrix.tar.gz

# Create a conda environment with baredSC:
wget https://raw.githubusercontent.com/lldelisle/baredSC/v1.0.0/baredSC_env.yml
conda env create -f baredSC_env.yml
conda activate baredSC

# Install anndata to be able to load h5ad:
conda install -y -c bioconda -c conda-forge anndata
conda deactivate
