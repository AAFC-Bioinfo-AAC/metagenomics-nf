#!/bin/bash -l
#SBATCH --job-name=metagenomics_nf
#SBATCH --output=stdout.out
#SBATCH --error=stderr.err
#SBATCH --partition=standard
#SBATCH --account=grdi_genarcc
#SBATCH --time=50:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --clusters gpsc7,gpsc8
#SBATCH --export=ALL

# This is a gold standard script that work well"

export PATH=/gpfs/fs7/aafc/pilot/aafc_sjsr/abcc/resources/miniconda3/bin/:$PATH
source /gpfs/fs7/aafc/pilot/aafc_sjsr/abcc/resources/miniconda3/etc/profile.d/conda.sh
conda activate nextflow

cd /gpfs/fs7/grdi/genarcc/wp2/brouardjs/metagenomic_nf

export NXF_OPTS="-Xms500M -Xmx2G"
export NXF_OFFLINE=true
export APPTAINER_TMPDIR=/fs/vnas_Haafc/qc/jsb000/tmp

nextflow run main.nf -profile gpsc -resume


