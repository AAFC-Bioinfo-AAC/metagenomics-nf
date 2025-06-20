## Setting up the .env file
The .env file is a text file essential to the pipeline operation. It contains various parameters necessary for the execution of the processes. These parameters are mainly paths to databases and slurm settings specific to your HPC. Make sure to fill in all parameters before running the analysis. A [template](./.env) is provided in the project directory.

- CONDA_SRC:  The path to the conda.sh file of your miniconda3 installation. Typically, it is `~/miniconda3/etc/profile.d/conda.sh`
- NXF_APPTAINER_CACHEDIR: Nextflow caches Apptainer images in the `apptainer` or `singularity` directory, in the pipeline work directory, by default. Set this to the same value as WORKDIR to be on the safe side.
- SLURM_ACCT: Your slurm account
- PARTITION: The partition for tasks that require less than 512GB
- PARTITION_HIGH: Set this parameter to a partition with more than 1TB of memory
- CLUSTERS: Specify the slurm cluster
- APPTAINER_IMGS: The location of your apptainer images
- WORKDIR: The location of Nextflow's `work` directory 
- GENOME: Basename of a reference genome for decontamination
- The remaining parameters are databases location of:  
CHOCOPHLAN_DB  
UNIREF_DB  
METAPHLAN_DB  
CHECKM2_DB  
GTDB_DB  
KAIJU_DB  
KRAKEN2_DB  
PHYLO_DB  

