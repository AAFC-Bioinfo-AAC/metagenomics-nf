# metagenomic_nf


## 1 - Introduction
Write an introduction here...

You can have a look at the [Workflow diagram](misc/flowchart.png).

## 2 - Dependencies


### 2.1 - Nextflow

This pipeline is build with the Nextflow language. If you are not familiar with this, you can read more [here](https://www.nextflow.io/). The [documentaton](https://www.nextflow.io/docs/latest/) is well written and updated regularly. A good approach is to install Nextflow itself in a conda environment.

### 2.2 - Conda environments

Nextflow and conda environments work well together. Almost all process of the workflow use specific conda environments.

However, creating multiple conda environments is tedious task. If you are working on the AAFC Biocluster, you can use pre-build conda environments that lie at this location :

```shell
/isilon/common/conda/brouardjs
```

Note that by default, these conda env will be used when running the pipeline with the **-profile biocluster** option.

Details about the recipes used to create these conda envs are detailed in the [conda_env.md](./conda_env.md) file.


### 2.3 - Databases

Several databases are required to perform all steps of the pipeline.

If you are running this pipeline on the AAFC Biocluster, you can take advantage of the pre-build databases whom location are specified in the Nextflow.config file.

You can read more about how databases were set-up in our [databases documentation](./databases.md).


## 3 - Prerequisites


### 3.1 - Preparation of a reference genome for decontamination


You will need to build an index of the host genome plus the phiX genome (Included in data/genome folder).

(Write basic instructions to dowload pig or cow genomes and cat them with the PhiX...)


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


When setting the rename parameters to 'yes', the rename workflow will rename the file id according to a map file that should be placed in the metadata folder.

The map file is a simple tsv file build with the raw sequences names. **Note that there are no header in tsv file!**


| id            |       Basenname of fastq file                    |
|---------------|--------------------------------------------------|
|C284-d21_R1	| NS.2055.003.IDT_i7_53---IDT_i5_53.C284-d21-WRC_R1|
|C284-d21_R2	| NS.2055.003.IDT_i7_53---IDT_i5_53.C284-d21-WRC_R2|
|C9009-d7_R1	| NS.2055.003.IDT_i7_55---IDT_i5_55.C9009-d7-WRC_R1|
|C9009-d7_R2	| NS.2055.003.IDT_i7_55---IDT_i5_55.C9009-d7-WRC_R2|
|C283-d28_R1	| NS.2055.003.IDT_i7_56---IDT_i5_56.C283-d28-WRC_R1|
|C283-d28_R2	| NS.2055.003.IDT_i7_56---IDT_i5_56.C283-d28-WRC_R2|



The following snippet will produce the desired map_file (Please adjust accordingly)


```shell
cd data/reads
for i in `ls *.fastq.gz`; do n=$(basename $i ".fastq.gz"); id=$(echo $i | cut -f 5 -d '.'); printf "$id\t$n\n"; done > ../../metadata/map_file.tsv
```


## 6 - Launch :rocket:

### 6.1 - AAFC Biocluster

```shell
# No resume
screen -S Run
conda activate nextflow
nextflow run main.nf -profile biocluster -resume ID 2>&1 | tee logfile_nextflow.txt

nextflow run main.nf -profile biocluster -resume d3bda63b-ed9d-4728-9b68-8171422cac65  2>&1 | tee logfile_nextflow_b.txt




# With resume and specifying the location of the work folder
nextflow run main.nf -profile biocluster -resume -w /isilon/projects/J-002888_GRDI-AMR2/work






```






### 6.2 - Waffles (NML cluster)

Don't use screen!!


export NXF_OPTS="-Xms500M -Xmx2G" 


```shell
conda activate nextflow-jsb
export NXF_OPTS="-Xms500M -Xmx2G" 
sbatch -D $PWD --export=ALL -J metagenomics_nf -c 2 --mem 4G -p NMLResearch -o $PWD/nextflow_log-%j.out --wrap="nextflow run main.nf -profile waffles -resume d370824e-10e5-464c-9628-4bb319978b55"
```


sbatch -D $PWD --export=ALL -J metagenomics_nf -c 2 --mem 4G -p NMLResearch -o $PWD/nextflow_log-%j.out --wrap="nextflow run main.nf  -profile waffles -ansi-log false -resume d370824e-10e5-464c-9628-4bb319978b55"




## 7 - Alternative paths


when looking at the whole pipeline, one can see that there are 3 main branches : main, humann and taxonomic profiling.

If you have your prepared_readsin hands, you can run easily run just one branch at a time by removing the undesired modules in the main.nf file. Here a screenshot of a sucessufully run of the kaiju branch!

<p align="center"><img src="misc/just_kaiju.png" alt="The kaiju branch"></p>



### In resume mode

You can use the resume command with the session ID to recover a specific execution. For example:

```shell
conda activate nextflow-jsb
# ceci a fonctionné avec la nouvelle branche orthodox!
sbatch -D $PWD --export=ALL -J metagenomics_nf -c 2 --mem 4G -p NMLResearch --wrap="nextflow run main.nf -profile waffles -resume 9025af6d-b3a6-4c63-bcf3-a98b6ee671d8 -with-report my_report"
```





