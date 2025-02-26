#!/bin/bash
#SBATCH --job-name=metagenomics-test
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --partition standard
#SBATCH --account=aafc_pilot
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
export PATH=/gpfs/fs7/aafc/pilot/aafc_sjsr/abcc/resources/miniconda3/bin/:$PATH
source /gpfs/fs7/aafc/pilot/aafc_sjsr/abcc/resources/miniconda3/etc/profile.d/conda.sh
conda activate nextflow
cd /gpfs/fs7/aafc/pilot/aafc_sjsr/abcc/metagenomics/metagenomic_nf
export NXF_OFFLINE=true
nextflow run main.nf -profile gpsc -resume