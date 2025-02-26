# Metagenomic_nf
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## ABOUT
This Nextflow workflow automates many metagenomics analyses from the quality filtering to the annotation of metagenomics assembled genomes (MAGs).
The pipeline includes several state-of-the-art programs in the field of metagenomics such as dRep, CheckM2, quast, phylophlan, DRAM, etc!

It is easy to dive into the code of this project as it includes three main files :

  * **nextflow.config** allows one to specify the profile and parameters of the analyses, the computational requirements of each task, etc.

  * **main.nf** describes the workflow logic.

  * **modules.nf** contains the Nextflow code for each bioinformatics program.

It has been proven to run seamlessly in at least two computing environments : the AAFC/AAC Biocluster and the PHAC/NML waffles cluster.

Boolean options allow the user to include or skip some components of the workflow : the Kaiju branch, the Kraken2/Bracken branch, the co-assembly branch.

The pipeline offers alternative entry points : it offers the possibility to start from prepared reads (reads that have been trimmed and decontaminated) or to specify already obtainded individual assemblies (with Megahit).

You can have a look at the [Workflow diagram](docs/misc/flowchart.png).

---

## TABLE OF CONTENTS
| **Section**                                | **Description**                                                                                           |
|--------------------------------------------|-----------------------------------------------------------------------------------------------------------|
| [ABOUT](#about)                            | A description of the Metagenomic_nf pipeline, introducing its purpose and process. |
| [OVERVIEW](#overview)                      | A detailed diagram of the pipeline's process to illustrate the workflow. |
| [DATA](#data)                              | Information about the input data used. |
| [PARAMETERS](#parameters)                  | A reference to the parameters used throughout the project. |
| [USAGE](#usage)                            | Instructions on setting up and running the pipeline. |
| &nbsp;&nbsp;&nbsp;&nbsp;[1. Dependencies](#1-dependencies) | Installation and configuration of necessary dependencies. |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[1.1 - Nextflow](#11---nextflow) | Introduction to Nextflow and installation instructions. |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[1.2 - Conda environments](#12---conda-environments) | Setting up conda environments for workflow execution. |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[1.3 - Databases](#13---databases) | Information on required databases and their setup. |
| &nbsp;&nbsp;&nbsp;&nbsp;[2. Prerequisites](#2-prerequisites) | Preliminary setup steps before running the pipeline. |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[2.1 - Preparation of a reference genome for decontamination](#21---preparation-of-a-reference-genome-for-decontamination) | Steps for preparing a Bowtie2 index of host and PhiX genome. |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[2.2 - Preparation of a map file](#22---preparation-of-a-map-file) | Instructions for generating a sample read-to-file mapping file. |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[2.3 - Preparation of a co-assembly map file](#23---preparation-of-a-co-assembly-map-file) | Creating a mapping file for co-assembly runs. |
| &nbsp;&nbsp;&nbsp;&nbsp;[3. Running the Pipeline](#3-running-the-pipeline) | Execution commands and configurations for different clusters. |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.1 - First launch](#31---first-launch) | Instructions for first-time pipeline execution. |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.1.1 - On AAFC AAC Biocluster](#311---on-aafc-aac-biocluster) | Running the pipeline on the AAFC/AAC Biocluster. |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.1.2 - On Waffles PHAC NML cluster](#312---on-waffles-phac-nml-cluster) | Running the pipeline on the Waffles cluster. |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.1.3 - On GPSC](#313---on-gpsc) | Running the pipeline on the GPSC cluster. |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.2 - Resume a run](#32---resume-a-run) | Resuming an interrupted or failed execution. |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.2.1 - On AAFC AAC Biocluster](#321---on-aafc-aac-biocluster) | Resume execution on AAFC/AAC Biocluster. |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.2.2 - On Waffles PHAC NML cluster](#322---on-waffles-phac-nml-cluster) | Resume execution on Waffles cluster. |
| &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[3.2.3 - On GPSC](#323---on-gpsc) | Resume execution on GPSC cluster. |
| [OUTPUT](#output)                          | A description of the output generated. |
| [CREDITS](#credits)                        | Acknowledgment of contributors, teams, and organizations that supported the project. |
| [CONTRIBUTION](#contribution)              | Guidelines for contributing to the repository, with a link to the `CONTRIBUTING.md` file. |
| [COPYRIGHT](#copyright)                    | Ownership details. |
| [LICENSE](#license)                        | Information about the license, including a link to the `LICENSE` file. |
| [PUBLICATIONS & ADDITIONAL RESOURCES](#publications-and-additional-resources) | Links to publications/articles/other resources related to the project. |
| [CITATION](#citation)                      | Instructions for citing the project, with references to the `CITATION.cff` and `CITATIONS.md` files. |

---

## OVERVIEW
<p align="center">
    <img src="./docs/misc/flowchart.png" alt="Metagenomics_nf flowchart." width="600" height="400">
</p>

---

## DATA

The pipeline processes metagenomic sequencing data. The different input data types include:

- **Raw sequencing reads** (paired-end, FASTQ format)  
  - These reads undergo quality control, trimming, and host contamination removal before downstream analyses.  
- **Prepared reads** (pre-processed reads after quality control)  
  - Users may provide already cleaned and decontaminated reads to bypass the initial preprocessing steps.  
- **Individual assemblies** (assembled contigs from tools like MEGAHIT)  
  - The workflow supports skipping the assembly step by accepting pre-assembled contigs.  
- **Reference genomes** (for host genome decontamination)  
  - Bowtie2 indexes of host genomes (e.g., pig, cow) are required for removing host DNA contamination.  
- **Database files** (for functional, taxonomic, and phylogenetic profiling)  
  - Several external databases are used, including:
    - **Kraken2 and Kaiju** (taxonomic classification)
    - **GTDB-Tk and PhyloPhlAn** (phylogenetic analysis)
    - **HUMAnN and DRAM** (functional profiling)
    - **CheckM2 and dRep** (MAGs quality control and dereplication)

---

## PARAMETERS

Pipeline execution is configured using the `nextflow.config` file and command-line parameters. Key configurable parameters include:

### **General Parameters**
- `--reads`  
  - Path to raw sequencing reads (`*_R1.fastq.gz`, `*_R2.fastq.gz`).
- `--prepared_reads`  
  - Path to pre-processed reads, if available.
- `--indiv_assemblies`  
  - Path to pre-assembled contigs (e.g., MEGAHIT output).
- `--genome`  
  - Path to the reference genome for host read filtering.
- `--genome_basename`  
  - Base name of the reference genome index used with Bowtie2.
- `--map_file`  
  - Sample metadata file mapping sample names to sequencing files.
- `--coassembly_file`  
  - Metadata file defining sample groups for co-assembly.

### **Workflow Execution Options**
- `--skip_kraken`  
  - Skip taxonomic classification using Kraken2 (default: false).
- `--skip_humann`  
  - Skip functional profiling with HUMAnN (default: false).
- `--skip_kaiju`  
  - Skip taxonomic classification using Kaiju (default: false).
- `--skip_coassembly`  
  - Skip co-assembly of metagenomic samples (default: false).
- `--use_prepared_reads`  
  - Use pre-filtered reads instead of raw reads (default: false).
- `--use_megahit_individual_assemblies`  
  - Use precomputed individual assemblies instead of running MEGAHIT (default: false).

### **Resource Allocation**
- `--cpus`  
  - Number of CPU cores allocated per process.
- `--memory`  
  - Amount of memory allocated per process.
- `--profile`  
  - Specifies the computing environment profile (`biocluster`, `waffles`, `gpsc`, or `local`).

---

## USAGE

### 1. Dependencies

#### 1.1 - Nextflow

This pipeline is built with the Nextflow language. If you are not familiar with this, you can read more [here](https://www.nextflow.io/). The [documentaton](https://www.nextflow.io/docs/latest/) is well written and updated regularly. A good approach is to install Nextflow itself in a conda environment.

#### 1.2 - Conda environments

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


#### 1.3 - Databases

Several databases are required to perform all steps of the pipeline.

***(todo: make a table with all databases and their versions)***

If you are running this pipeline on the AAFC Biocluster, the PHAC/NML waffles cluster, or the GPSC cluster, you can take advantage of the pre-built databases whose location are specified in the Nextflow.config file.

You can read more about how databases were set-up in our [databases documentation](./databases.md).


### 2. Prerequisites

#### 2.1 - Preparation of a reference genome for decontamination


You will need to build a Bowtie2 index of the host genome plus the phiX genome (Included in data/genome folder).

Note that if you are working on the PHAC/NML waffles cluster, pre-built indexes are available at this location :

```shell
/Drives/O/GRDI-AMR2/share/genomes
```

##### 2.1.1 - Example with the pig genome

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

#### 2.2 - Preparation of a map file

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

#### 2.3 - Preparation of a co-assembly map file

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

### 3. Running the Pipeline:

#### 3.1 - First launch

##### 3.1.1 - On AAFC AAC Biocluster

You simply use the biocluster profile:

```shell
# Warning!  No resume!
screen -S Run
conda activate nextflow
nextflow run main.nf -profile biocluster | tee logfile_nextflow.txt
```

You can detach of your screen using ctlr + a + d.

##### 3.1.2 - On Waffles PHAC NML cluster

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
##### 3.1.3 - On GPSC
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

#### 3.2 - Resume a run

##### 3.2.1 - On AAFC AAC Biocluster

You can use the resume command with the session ID to recover a specific execution. For example:

```shell
# Obtain the desired run id
nextflow log
nextflow run main.nf -profile biocluster -resume d3bda63b-ed9d-4728-9b68-8171422cac65  2>&1 | tee logfile_nextflow.txt
```

##### 3.2.2 - On Waffles PHAC NML cluster

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
##### 3.2.3 - On GPSC
#Insert info here

---

## OUTPUT
Pipeline outputs are directed to the folder specified in the results parameter from the profile used in `nextflow.config`. By default this is the output directory.

---

## FUTURE DIRECTIONS

Our intention is to include other modules related to the detection of AMR genes and plasmidic sequences. Our list includes ABRIcate with MEGAres v3, PlasForest, Mob-Suite, geNomad, mobileOG-db and Sarand.

---

## KNOWN ISSUES

Feel free to use the Issues section to report any problems you may encounter with this pipeline.

### Storage burden issue

When used with a large number of large samples, the work folder can grow up to several terabytes in spite of the fact that the workflows use diverse strategies to mitigate the number and size of temporary files.

A good strategy would be to not run all of the different branches of the workflow at the beginning. For example, you can begin by skipping the taxonomic and co-assembly branches, trying to get **all your prepared reads** and **all your individual assemblies**. Upon completion, copy them from the results folder (cp -L) to a secure place. Then remove the large work folder and rerun the workflow by specifying the location of your prepared reads and individual assemblies in the config file.

### Co-assembly issue

When run with a large number of metagenomic samples (50-100), the co-assembly can be very complex and can take **weeks** to complete! To overcome this issue you have the option to perform several co-assemblies each containing a reduced number of samples (~2-10 samples).

---

##  CREDITS

The metagenomic_nf workflow was written in the Nextflow language by Jean-Simon Brouard (AAFC/AAC Sherbrooke RDC). The main components of this workflow come from the work of Devin Holman (AAFC/AAC Lacombe RDC) whereas the original scripts were written in Bash by Arun Kommadath (AAFC/AAC Lacombe). Sara Ricci, from the team of Renee Petri (AAFC/AAC Sherbrooke RDC) has also contributed to adapt this workflow for being used with cow samples.

---

## CONTRIBUTION
If you would like to contribute to this pipeline, please see the [contributing guidelines](https://gccode.ssc-spc.gc.ca/abcc_rcba/documentation-guide/-/blob/main/Guides/how_to_contribute.md?ref_type=heads#contributing-to-someone-elses-pipeline).

---

## COPYRIGHT
Government of Canada, Agriculture & Agri-Food Canada

---

## LICENSE 
This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## PUBLICATIONS AND ADDITIONAL RESOURCES

* [Conda Environment Setup](docs/conda_env.md)
* [Database Information](docs/databases.md)
* [Workflow Document](docs/Metagenomic_workflow_March_9_2023_Ricci.docx)

The Metagenomic_nf workflow integrates various state-of-the-art tools for metagenomics analysis. Below are key references and resources for the tools used in the pipeline:

### **1. General Nextflow Documentation**
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [Nextflow GitHub Repository](https://github.com/nextflow-io/nextflow)

### **2. Quality Control**
- FastQC: [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- MultiQC: [https://multiqc.info/](https://multiqc.info/)

### **3. Taxonomic Classification**
- Kraken2: [https://ccb.jhu.edu/software/kraken2/](https://ccb.jhu.edu/software/kraken2/)
- Bracken: [https://ccb.jhu.edu/software/bracken/](https://ccb.jhu.edu/software/bracken/)
- Kaiju: [https://github.com/bioinformatics-centre/kaiju](https://github.com/bioinformatics-centre/kaiju)

### **4. Assembly & Binning**
- MEGAHIT: [https://github.com/voutcn/megahit](https://github.com/voutcn/megahit)
- MetaBAT2: [https://bitbucket.org/berkeleylab/metabat/src/master/](https://bitbucket.org/berkeleylab/metabat/src/master/)
- QUAST: [http://bioinf.spbau.ru/quast](http://bioinf.spbau.ru/quast)

### **5. Functional Profiling**
- HUMAnN: [https://huttenhower.sph.harvard.edu/humann/](https://huttenhower.sph.harvard.edu/humann/)
- DRAM: [https://github.com/shafferm/DRAM](https://github.com/shafferm/DRAM)

### **6. Phylogenetic Analysis**
- GTDB-Tk: [https://github.com/Ecogenomics/GTDBTk](https://github.com/Ecogenomics/GTDBTk)
- PhyloPhlAn: [https://github.com/biobakery/phylophlan](https://github.com/biobakery/phylophlan)

### **7. Dereplication & MAG Refinement**
- CheckM2: [https://github.com/chklovski/CheckM2](https://github.com/chklovski/CheckM2)
- dRep: [https://github.com/MrOlm/drep](https://github.com/MrOlm/drep)

For additional inquiries or troubleshooting, consult the **Issues** section of this repository.

---

## CITATION
If you use this pipeline for your analysis, please cite it using the information found in [CITATION.cff](CITATION.cff).

An extensive list of references for the tools used by the pipeline can be found in the [CITATIONS.md](LINK) file.