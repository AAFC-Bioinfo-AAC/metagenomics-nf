#!/bin/bash

source .env
source ${CONDA_SRC}
export NXF_OPTS="-Xms500M -Xmx2G"
export NXF_OFFLINE=true
export NXF_APPTAINER_CACHEDIR=${NXF_APPTAINER_CACHEDIR}

conda activate nextflow

sbatch -D $PWD \
  --account=${SLURM_ACCT} \
  --partition=${PARTITION} \
  --clusters=${CLUSTERS} \
  --job-name=metagenomic_nf \
  --output=results/logs/stdout.out \
  --error=results/logs/stderr.err \
  --time=48:00:00 \
  --ntasks=1 \
  --nodes=1 \
  --cpus-per-task=2 \
  --export=ALL \
  --wrap="nextflow -log ./results/logs/nextflow.log run ./main.nf -profile hpc ./nextflow.config -resume"

