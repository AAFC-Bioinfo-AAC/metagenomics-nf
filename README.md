# metagenomic_nf


## 1 - Introduction

This Nextflow workflow allows the realization of the most common metagenomics analyses from the quality filtering to the annotation of metagenomics assembled genomes (MAGs).

The pipeline includes several state-of-the-art programs in the field of metagenomics such as dRep, CheckM2, quast, phylophlan, DRAM, etc!

### Easy to set-up/adapt

It is easy to dive into the code of this projects as it includes three main files :

  * **nextflow.config** allows one to specify the profile and parameters of the analyses, the computational requirements of each task, etc.

  * **main.nf** describes the workflow logic.

  * **modules.nf** contains the Nextflow code for each bioinformatics programs.

### Versatile

It has proven to run seemlesly in a least to 2 computing environments : the AAFC/AAC Biocluster and the PHAC/NML waffles cluster.

Boolean options allows the users to include or skip some components of the workflow : the Kaiju branch, the Kraken2/Bracken branch, the co-assembly branch.

The pipeline offers alternative entry points : it offers the possibility to start from prepared reads (reads that have been trimmed and decontaminated) or to specify already obtainded individual assemblies (With Megahit).


You can have a look at the [Workflow diagram](misc/flowchart.png).

## 2 - Dependencies


### 2.1 - Nextflow

This pipeline is build with the Nextflow language. If you are not familiar with this, you can read more [here](https://www.nextflow.io/). The [documentaton](https://www.nextflow.io/docs/latest/) is well written and updated regularly. A good approach is to install Nextflow itself in a conda environment.

### 2.2 - Conda environments

Nextflow and conda environments work well together. Almost all process of the workflow use specific conda environments.

However, creating multiple conda environments is a tedious task. If you are working on the AAFC Biocluster, you can use the pre-build conda environments that lie at this location :

```shell
/isilon/common/conda/brouardjs
```

By default, these conda env will be used when running the pipeline with the **-profile biocluster** option.

Details about the recipes used to create these conda envs are detailed in the [conda_env.md](./conda_env.md) file.


### 2.3 - Databases

Several databases are required to perform all steps of the pipeline.

***(todo: make a table with all databases and their versions)***

If you are running this pipeline on the AAFC Biocluster, you can take advantage of the pre-build databases whom location are specified in the Nextflow.config file.

You can read more about how databases were set-up in our [databases documentation](./databases.md).


## 3 - Prerequisites


### 3.1 - Preparation of a reference genome for decontamination


You will need to build a Bowtie2 index of the host genome plus the phiX genome (Included in data/genome folder).

If you are using cow or pig samples, you can throw me an email and I could share my pre-build indexes.


***(to do: write basic instructions to dowload pig or cow genomes and cat them with the PhiX...)***


### 3.1.1 - Example with the pig genome



**Waffles cluster (NML) (slurm)**

```shell
conda activate bowtie2
sbatch \
 -D $PWD \
 --output $PWD/bowtie2-$j.out \
 --export=ALL \
 -J bowtie2-build \
 -c 8 \
 -p NMLResearch \
 --mem 64G \
 --wrap="bowtie2-build Pig_PhiX_genomes.fna pig/pig"
```


### 3.1.2 - Example with the cow genome


...



## 3.2 - Preparation of a map file


When setting the rename parameters to 'yes', the rename sub-workflow will rename the file id according to a map file that should be placed in the metadata folder.

The map file is a simple tsv file build with the raw sequences names. **Note that there are no header in tsv file!**


| id            |       Basenname of fastq file                    |
|---------------|--------------------------------------------------|
|C284-d21_R1	| NS.2055.003.IDT_i7_53---IDT_i5_53.C284-d21-WRC_R1|
|C284-d21_R2	| NS.2055.003.IDT_i7_53---IDT_i5_53.C284-d21-WRC_R2|
|C9009-d7_R1	| NS.2055.003.IDT_i7_55---IDT_i5_55.C9009-d7-WRC_R1|
|C9009-d7_R2	| NS.2055.003.IDT_i7_55---IDT_i5_55.C9009-d7-WRC_R2|
|C283-d28_R1	| NS.2055.003.IDT_i7_56---IDT_i5_56.C283-d28-WRC_R1|
|C283-d28_R2	| NS.2055.003.IDT_i7_56---IDT_i5_56.C283-d28-WRC_R2|



The following snippet will produce the desired map_file (Please adjust accordingly) using your raw input files!


```shell
cd data/reads
for i in `ls *.fastq.gz`; do n=$(basename $i ".fastq.gz"); id=$(echo $i | cut -f 5 -d '.'); printf "$id\t$n\n"; done > ../../metadata/map_file.tsv
```


## 4 - Launch :rocket:

### 4.1 - Launch on the AAFC/AAC Biocluster

```shell
# No resume
screen -S Run
conda activate nextflow

# No resume
nextflow run main.nf -profile biocluster 2>&1 | tee logfile_nextflow.txt

# Specifying the location of the work folder
nextflow run main.nf -profile biocluster -w /isilon/projects/J-002888_GRDI-AMR2/work
```

### 4.2 - Resume a run on the AAFC Biocluster

You can use the resume command with the session ID to recover a specific execution. For example:

```shell
# Obtain the desired run id
nextflow log

# Specif a specific run to resume
nextflow run main.nf -profile biocluster -resume d3bda63b-ed9d-4728-9b68-8171422cac65  2>&1 | tee logfile_nextflow.txt
```


### 4.3 - Launch on Waffles (PHAC/NML cluster)

Don't use screen!!

```shell
conda activate nextflow-jsb
export NXF_OPTS="-Xms500M -Xmx2G" 
sbatch -D $PWD --export=ALL -J metagenomics_nf -c 2 --mem 4G -p NMLResearch -o $PWD/nextflow_log-%j.out --wrap="nextflow run main.nf -profile waffles"
```


### 4.4 - Resume a run on Waffles (PHAC/NML cluster)

```shell
conda activate nextflow-jsb
export NXF_OPTS="-Xms500M -Xmx2G" 

# Obtain the desired run id
nextflow log

sbatch -D $PWD --export=ALL -J metagenomics_nf -c 2 --mem 4G -p NMLResearch -o $PWD/nextflow_log-%j.out --wrap="nextflow run main.nf -profile waffles -resume 94de4004-69dc-4a11-9cef-c936e89974a3"
```

## 5 - Future directions

Our intention is to include others modules related to the detection of AMR genes, plasmidic sequences. Our list includes ABRIcate with MEGAres v3, AMRplusplus, PlasForest, Mob-Suite, geNomad, mobileOG-db and Sarand.


## 6 - Issues

Feel free to use the Issue section to report any problems you may encounter with this pipeline.


### 6.1 - Storage burden issue

When used with a large number of large samples, the work folder can growth up to several terabytes in spite that the workflows uses diverse strategies to mitigate the number and size of temporary files.

A good strategy would be to not run all different branches of the workflow at the beginning. For example, you can begin by skipping the taxonomic and co-assembly branches, trying to get **all your prepared reads** and **all your individual assemblies**. Upon completion, copy them from the results foder (cp -L) to a secure place. Then remove the large work folder and rerun the workflow by specifying the location of your prepared reads and individual assemblies in the config file.

### 6.2 - Co-assembly issue

When run with a large number of metagenomic samples (50-100), the co-assembly can be very complex and can takes **weeks** to complete! The entire pipeline can be run by skipping this step. It's up to you to decide if it's worth it.



##  7 - Credits

Metagenomic_nf workflow was written in the Nextflow language by Jean-Simon Brouard (AAFC/AAC Sherbrooke RDC). The main components of this workflow comes from the work of Devin Holman (AAFC/AAC Lacombe RDC) whereas the original scripts were written in Bash by Arun Kommadath (AAFC/AAC Lacombe). Sara Ricci, from the team of Renee Petri (AAFC/AAC Sherbrooke RDC) has also contributed to adapt this workflow for being used with cow samples.

