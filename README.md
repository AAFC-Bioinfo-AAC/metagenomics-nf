# metagenomic_nf


## 1 - Introduction

This Nextflow workflow automates many metagenomics analyses from the quality filtering to the annotation of metagenomics assembled genomes (MAGs).

The pipeline includes several state-of-the-art programs in the field of metagenomics such as dRep, CheckM2, quast, phylophlan, DRAM, etc!

### Easy to set-up/adapt

It is easy to dive into the code of this project as it includes three main files :

  * **nextflow.config** allows one to specify the profile and parameters of the analyses, the computational requirements of each task, etc.

  * **main.nf** describes the workflow logic.

  * **modules.nf** contains the Nextflow code for each bioinformatics programs.

### Versatile

It has been proven to run seamlessly in at least two computing environments : the AAFC/AAC Biocluster and the PHAC/NML waffles cluster.

Boolean options allow the user to include or skip some components of the workflow : the Kaiju branch, the Kraken2/Bracken branch, the co-assembly branch.

The pipeline offers alternative entry points : it offers the possibility to start from prepared reads (reads that have been trimmed and decontaminated) or to specify already obtainded individual assemblies (with Megahit).


You can have a look at the [Workflow diagram](misc/flowchart.png).

## 2 - Dependencies


### 2.1 - Nextflow

This pipeline is built with the Nextflow language. If you are not familiar with this, you can read more [here](https://www.nextflow.io/). The [documentaton](https://www.nextflow.io/docs/latest/) is well written and updated regularly. A good approach is to install Nextflow itself in a conda environment.

### 2.2 - Conda environments

Nextflow and conda environments work well together. Almost all processes of the workflow use specific conda environments.

However, creating multiple conda environments is a tedious task.

If you are working on the AAFC Biocluster, you can use the pre-built conda environments that lie at this location :

```shell
/isilon/common/conda/brouardjs
```

If you are working on the PHAC/NML waffles cluster, you can use the pre-built conda environments that lie at this location :

```shell
/Drives/O/GRDI-AMR2/share/conda
```

If you are working on the GPSC cluster, you can use the pre-built conda environments that lie at this location :

```shell
/gpfs/fs7/aafc/pilot/aafc_sjsr/abcc/resources/env/metagenomics_nf
```

By default, these conda envs will be used when running the pipeline with the **-profile biocluster**, **-profile waffles**, or **--profile gpsc** options.

Details about the recipes used to create these conda envs are detailed in the [conda_env.md](./conda_env.md) file.


### 2.3 - Databases

Several databases are required to perform all steps of the pipeline.

***(todo: make a table with all databases and their versions)***

If you are running this pipeline on the AAFC Biocluster, the PHAC/NML waffles cluster, or the GPSC cluster, you can take advantage of the pre-built databases whose location are specified in the Nextflow.config file.

You can read more about how databases were set-up in our [databases documentation](./databases.md).


## 3 - Prerequisites


### 3.1 - Preparation of a reference genome for decontamination


You will need to build a Bowtie2 index of the host genome plus the phiX genome (Included in data/genome folder).

Note that if you are working on the PHAC/NML waffles cluster, pre-built indexes are available at this location :

```shell
/Drives/O/GRDI-AMR2/share/genomes
```


#### 3.1.1 - Example with the pig genome


You can download the latest pig reference genome here:
```shell
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/025/GCF_000003025.6_Sscrofa11.1/GCF_000003025.6_Sscrofa11.1_genomic.fna.gz
```

Then unzip the file:

```shell
gunzip GCF_000003025.6_Sscrofa11.1_genomic.fna.gz
```

Combine the pig and PhiX genomes:

```shell
cat data/genomes/phiX.fa GCF_000003025.6_Sscrofa11.1_genomic.fna > Pig_PhiX_genomes.fna
```

Finally, build the Bowtie2 index:  
- On the Waffle cluster:
```shell
prefix='/Drives/O/GRDI-AMR2/share/conda/envs'
conda activate $prefix/bowtie2_2.5.2
mkdir pig
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
- On the GPSC:  
```shell
prefix='/gpfs/fs7/aafc/pilot/aafc_sjsr/abcc/resources/env/metagenomics_nf'
conda activate $prefix/bowtie2
mkdir pig
sbatch \
 -D $PWD \
 --output $PWD/bowtie2-$j.out \
 --export=ALL \
 -J bowtie2-build \
 -c 8 \
 -p standard \
 --account=aafc_pilot \
 -t 300 \ 
 --wrap="bowtie2-build Pig_PhiX_genomes.fna pig/pig"
 ```


### 3.2 - Preparation of a map file

The map file allows you to give a relevant sample name to the corresponding metagenomics reads.

The map file is a .tsv file with at least 2 mandatory columns : sample_read and file_name. The idea is to have an expression corresponding to the desired sampleID plus _R1 or 
_R2 in the first column and the basename of the corresponding fastq files in the second column. Here is an example :


| sample_read | file_name	|
|-------------|-----------------------------------------------|
|G2-E1-13_R1	|NS.2066.001.IDT_i7_13---IDT_i5_13.G2-E1-13_R1|
|G2-E1-13_R2	|NS.2066.001.IDT_i7_13---IDT_i5_13.G2-E1-13_R2|	
|G2-E1-38_R1	|NS.2066.001.IDT_i7_14---IDT_i5_14.G2-E1-38_R1|
|G2-E1-38_R2	|NS.2066.001.IDT_i7_14---IDT_i5_14.G2-E1-38_R2|	
|G2-E1-63_R1	|NS.2066.001.IDT_i7_15---IDT_i5_15.G2-E1-63_R1|	
|G2-E1-63_R2	|NS.2066.001.IDT_i7_15---IDT_i5_15.G2-E1-63_R2|
|G2-E1-28_R1	|NS.2066.001.IDT_i7_51---IDT_i5_51.G2-E1-28_R1|	
|G2-E1-28_R2	|NS.2066.001.IDT_i7_51---IDT_i5_51.G2-E1-28_R2|	
|G2-E1-16_R1	|NS.2066.001.IDT_i7_39---IDT_i5_39.G2-E1-16_R1|	
|G2-E1-16_R2	|NS.2066.001.IDT_i7_39---IDT_i5_39.G2-E1-16_R2|


When setting the rename parameters to 'yes', the renamed sub-workflow will rename the reads according to the values present in the sample_read column of the map file.

A good starting point to produce a compliant map_file for this workflow is to use the following snippet with your raw input files!
Depending on how your raw reads are named, you may have to adjust the cut command however.


```shell
cd data/reads
printf "sample_read\tfile_name\n" >  ../../metadata/map_file.tsv
for i in `ls *.fastq.gz`; do n=$(basename $i ".fastq.gz"); id=$(echo $i | cut -f 5 -d '.'); printf "$id\t$n\n"; done >> ../../metadata/map_file.tsv
```



### 3.3 - Preparation of a co-assembly map file

If you want to perform co-assemblies, you will have also to prepare **a distinct file** specifying the desired groups for co-assemblies, e.g. :



| sample_read	| file_name	| project	| sample_type	| coassembly_group |
|-------------|-----------|---------|-------------|-------------------|
|G2-E1-13_R1	|NS.2066.001.IDT_i7_13---IDT_i5_13.G2-E1-13_R1	|StressCuZnII	|feces	  |stress_feces|
|G2-E1-13_R2	|NS.2066.001.IDT_i7_13---IDT_i5_13.G2-E1-13_R2	|StressCuZnII	|feces	  |stress_feces|
|G2-E1-38_R1	|NS.2066.001.IDT_i7_14---IDT_i5_14.G2-E1-38_R1	|LLQ	        |feces	  |llq_feces_herd5|
|G2-E1-38_R2	|NS.2066.001.IDT_i7_14---IDT_i5_14.G2-E1-38_R2	|LLQ	        |feces	  |llq_feces_herd5|
|G2-E1-63_R1	|NS.2066.001.IDT_i7_15---IDT_i5_15.G2-E1-63_R1	|PH	          |sediment	|ph_sediment|
|G2-E1-63_R2	|NS.2066.001.IDT_i7_15---IDT_i5_15.G2-E1-63_R2	|PH	          |sediment	|ph_sediment|
|G2-E1-28_R1	|NS.2066.001.IDT_i7_51---IDT_i5_51.G2-E1-28_R1	|LLQ	        |feces	  |llq_feces_herd5|
|G2-E1-28_R2	|NS.2066.001.IDT_i7_51---IDT_i5_51.G2-E1-28_R2	|LLQ	        |feces	  |llq_feces_herd5|
|G2-E1-16_R1	|NS.2066.001.IDT_i7_39---IDT_i5_39.G2-E1-16_R1	|StressCuZnII	|feces	  |stress_feces|
|G2-E1-16_R2	|NS.2066.001.IDT_i7_39---IDT_i5_39.G2-E1-16_R2	|StressCuZnII	|feces	  |stress_feces|


The co-assembly map file has 3 mandatory columns : sample_read, file_name and coassembly_group.


## 4 - Launch :rocket:

### 4.1 - First launch

#### 4.1.1  On AAFC/AAC Biocluster

You simply use the biocluster profile:

```shell
# Warning!  No resume!
screen -S Run
conda activate nextflow
nextflow run main.nf -profile biocluster | tee logfile_nextflow.txt
```

You can detach of your screen using ctlr + a + d.

#### 4.1.2  On Waffles (PHAC/NML cluster)

You simply use the waffles profile:

Don't use screen!

```shell
conda activate your-nextflow-env
export NXF_OPTS="-Xms500M -Xmx2G" 
sbatch -D $PWD \
       --export=ALL \
       -J metagenomics_nf \
       -c 2 \
       --mem 4G \
       -p NMLResearch \
       -o $PWD/nextflow_log-%j.out \
       --wrap="nextflow run main.nf -profile waffles -resume"
```
#### 4.1.3 On GPSC
You simply use the GPSC profile:

```shell
conda activate nextflow
export NXF_OPTS="-Xms500M -Xmx2G" 
sbatch -D $PWD \
       --export=ALL \
       -J metagenomics_nf \
       -c 2 \
       --mem 4G \
       -p standard \
       --account=aafc_pilot \
       -t 300 \
       -o $PWD/nextflow_log-%j.out \
       --wrap="nextflow run main.nf -profile gpsc"
```

### 4.2 - Resume a run

#### 4.2.1  On AAFC/AAC Biocluster

You can use the resume command with the session ID to recover a specific execution. For example:

```shell
# Obtain the desired run id
nextflow log
nextflow run main.nf -profile biocluster -resume d3bda63b-ed9d-4728-9b68-8171422cac65  2>&1 | tee logfile_nextflow.txt
```


#### 4.2.2 On Waffles (PHAC/NML cluster)

```shell
conda activate your-nextflow-env
export NXF_OPTS="-Xms500M -Xmx2G" 

# Obtain the desired run id
nextflow log

sbatch -D $PWD \
       --export=ALL \
      -J metagenomics_nf \
      -c 2 \
      --mem 4G \
      -p NMLResearch \
      -o $PWD/nextflow_log-%j.out \
      --wrap="nextflow run main.nf -profile waffles -resume 94de4004-69dc-4a11-9cef-c936e89974a3"
```
#### 4.2.2 On GPSC
#Insert info here

## 5 - Future directions

Our intention is to include other modules related to the detection of AMR genes and plasmidic sequences. Our list includes ABRIcate with MEGAres v3, PlasForest, Mob-Suite, geNomad, mobileOG-db and Sarand.


## 6 - Issues

Feel free to use the Issues section to report any problems you may encounter with this pipeline.


### 6.1 - Storage burden issue

When used with a large number of large samples, the work folder can grow up to several terabytes in spite of the fact that the workflows use diverse strategies to mitigate the number and size of temporary files.

A good strategy would be to not run all of the different branches of the workflow at the beginning. For example, you can begin by skipping the taxonomic and co-assembly branches, trying to get **all your prepared reads** and **all your individual assemblies**. Upon completion, copy them from the results folder (cp -L) to a secure place. Then remove the large work folder and rerun the workflow by specifying the location of your prepared reads and individual assemblies in the config file.

### 6.2 - Co-assembly issue

When run with a large number of metagenomic samples (50-100), the co-assembly can be very complex and can take **weeks** to complete! To overcome this issue you have the option to perform several co-assemblies each containing a reduced number of samples (~2-10 samples).


##  7 - Credits

The metagenomic_nf workflow was written in the Nextflow language by Jean-Simon Brouard (AAFC/AAC Sherbrooke RDC). The main components of this workflow come from the work of Devin Holman (AAFC/AAC Lacombe RDC) whereas the original scripts were written in Bash by Arun Kommadath (AAFC/AAC Lacombe). Sara Ricci, from the team of Renee Petri (AAFC/AAC Sherbrooke RDC) has also contributed to adapt this workflow for being used with cow samples.

