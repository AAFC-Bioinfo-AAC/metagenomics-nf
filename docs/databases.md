# Metagenomic_nf databases preparation

This guide provides step-by-step instructions for downloading, configuring, and verifying key databases used in this metagenomic workflow. Each section corresponds to a specific tool or package, including Kaiju, CheckM2, HUMAnN3, and others. These instructions are tailored for reproducible research, assuming the use of [Conda environments](./conda_env.md) and Linux-based systems. Be sure to adjust paths and environment variables as needed for your system.

## 1 - Kaiju
Kaiju is used for taxonomic classification of metagenomic sequences based on translated protein-level comparisons. Below are the commands to build the Kaiju database using the `nr_euk` index. Full Instructions can be found [here](https://github.com/bioinformatics-centre/kaiju).

Building the kaiju database requires a lot of memory (~400GB). If you're using an HPC system, the process may need to be split into two steps if the interactive node available has not enough memory.

Start by downloading and extracting the required files:
```shell
conda activate kaiju
kaiju-makedb -s nr_euk -t 16
```
If the compute node has enough memory, the command will complete successfully, and the database will be fully built. Allow ~**two days** of calculation time to complete the operation.

But if the compute node has not enough memory, you could get an output like this one:

```bash
Downloading taxdump.tar.gz
taxdump.tar.gz              100%[=========================================>]  66.32M  87.3MB/s    in 0.8s    
2025-04-23 12:39:26 URL:https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz [69537442/69537442] -> "taxdump.tar.gz" [1]
Extracting taxdump.tar.gz
Downloading nr.gz
nr.gz                       100%[=========================================>] 186.28G   105MB/s    in 24m 24s 
2025-04-23 13:03:58 URL:https://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz [200016146013/200016146013] -> "nr_euk/nr.gz" [1]
Downloading prot.accession2taxid.gz
prot.accession2taxid.gz     100%[=========================================>]   9.71G   147MB/s    in 76s     
2025-04-23 13:05:14 URL:https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz [10425285857/10425285857] -> "nr_euk/prot.accession2taxid.gz" [1]
Converting NR file to Kaiju database
13:05:14 Reading taxonomic tree from file nodes.dmp
13:05:16 Reading file merged.dmp
13:05:16 Reading taxa from file /home/user/miniconda3/envs/kaiju/bin/kaiju-taxonlistEuk.tsv
13:05:18 Reading accession to taxon id map from file nr_euk/prot.accession2taxid.gz
Killed
(kaiju) 
```
In this case, the downloads completed successfully, but the process was killed due to insufficient memory during the indexing step.

Assuming all required files were downloaded, you can finish the database build from those files by submitting a SLURM job on a partition with enough memory and with `--no-download` as an additional argument:
```bash
conda activate kaiju
sbatch \
 -D $PWD \
 --output $PWD/kaiju-%j.out \
 --export=ALL \
 -J kaiju_db_build \
 -c 16 \
 -p your_partition \
 --account=your_account \
 -t 4320 \
 --wrap="kaiju-makedb -s nr_euk -t 16 --no-download"
```
Kaiju only needs the files `nr_euk/kaiju_db_nr_euk.fmi`, `nodes.dmp` and `names.dmp` and the fmi file should be >130G. The remaining files can be deleted.


## 2 - Checkm2
[CheckM2](https://github.com/chklovski/CheckM2?tab=readme-ov-file#installation) is used for assessing the quality of microbial genomes. 

Create the database with:
```bash
conda activate checkm2
checkm2 database --download --path /absolute/path/to/database
```
After completion of the above command, a file named `CheckM2_database/uniref100.KO.1.dmnd` of ~3.1GB will appear in the current directory. Finally, CHECKM2_DB in the `.env` file must point to `/absolute/path/to/CheckM2_database/uniref100.KO.1.dmnd`  

## 3 - Humann3

### 3.1 - Download databases with the built in humann_databases command

HUMAnN3 requires several databases to function properly: pangenome, protein, and annotation databases. These should be downloaded using the `humann_databases` command. Make sure your `biobakery3` environment is active, and define a standard installation location before proceeding. 

#### 3.1.1 - Pangenome database

To download or upgrade the ChocoPhlAn pangenome database:

```bash
conda activate biobakery3
INSTALL_LOCATION=HUMAnN_db/
#To download or upgrade your pangenome database: 
humann_databases --download chocophlan full $INSTALL_LOCATION --update-config yes
```

#### 3.1.2 - Protein database

To download or update the UniRef90 protein database:
```bash
conda activate biobakery3
INSTALL_LOCATION=HUMAnN_db/
#To download or upgrade your protein database: 
humann_databases --download uniref uniref90_diamond $INSTALL_LOCATION --update-config yes
```

#### 3.1.3 - Annotation database

To download the utility mapping database used for annotation:
```bash
conda activate biobakery3
INSTALL_LOCATION=HUMAnN_db/
#To download or upgrade your ut database: 
humann_databases --download utility_mapping full $INSTALL_LOCATION --update-config yes
```


### 3.1.4 - Verify your installation
To ensure the databases were properly registered, you can check the location of these databases with :

```bash
conda activate biobakery3
humann_config --print
```
Output expected:
```
HUMAnN Configuration ( Section : Name = Value )
database_folders : nucleotide = /path/to/HUMAnN_db/chocophlan
database_folders : protein = /path/to/HUMAnN_db/uniref
database_folders : utility_mapping = /path/to/HUMAnN_db/utility_mapping
run_modes : resume = False
run_modes : verbose = False
run_modes : bypass_prescreen = False
run_modes : bypass_nucleotide_index = False
run_modes : bypass_nucleotide_search = False
run_modes : bypass_translated_search = False
run_modes : threads = 1
alignment_settings : evalue_threshold = 1.0
alignment_settings : prescreen_threshold = 0.01
alignment_settings : translated_subject_coverage_threshold = 50.0
alignment_settings : translated_query_coverage_threshold = 90.0
alignment_settings : nucleotide_subject_coverage_threshold = 50.0
alignment_settings : nucleotide_query_coverage_threshold = 90.0
output_format : output_max_decimals = 10
output_format : remove_stratified_output = False
output_format : remove_column_description_output = False
```


## 4 - Metaphlan

MetaPhlAn is used for profiling the composition of microbial communities based on clade-specific marker genes. However, due to version compatibility issues between MetaPhlAn and HUMAnN, special attention must be paid to which database version is installed.

### 4.1 - Fix an annoying bug that lead MetaPhlAn from trying to download the latest database!

If you’re using HUMAnN 3.6 with MetaPhlAn 4.0.3, you may encounter issues with automatic database downloads. MetaPhlAn may attempt to fetch a newer, incompatible version of the database. Here's an example of the error message:

```bash
  Downloading MetaPhlAn database
  Please note due to the size this might take a few minutes

  File "$METAPHLAN_DB" already present!

  Downloading http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212.md5
  Downloading file of size: 0.00 MB
  0.01 MB 11702.86 %  74.14 MB/sec  0 min -0 sec
  MD5 checksums do not correspond! If this happens again, you should remove the database files and rerun MetaPhlAn so they are re-downloaded
```
This is caused by MetaPhlAn trying to use the latest version by default. HUMAnN 3.6 only supports the _mpa_vJan21_CHOCOPhlAnSGB_202103_ version at this time. You can specify this index manually with the `--metaphlan-options` flag.

> More context and discussion:  
> [https://forum.biobakery.org/t/humann3-compatibility-problem-with-metaphlan4/4894/2](https://forum.biobakery.org/t/humann3-compatibility-problem-with-metaphlan4/4894/2)  
> [https://forum.biobakery.org/t/humann3-compatibility-problem-with-metaphlan4/4894/5](https://forum.biobakery.org/t/humann3-compatibility-problem-with-metaphlan4/4894/5)

### 4.2 - Solution: install the **latest compatible db** with Metaphlan3.6
To install the compatible version of the MetaPhlAn database at a custom location, run:

```bash
conda activate biobakery3
INSTALL_LOCATION=HUMAnN_db/
metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db $INSTALL_LOCATION/metaphlan/
```

If you're working on a high-performance computing cluster, use command like this:

```bash
conda activate biobakery3
INSTALL_LOCATION=HUMAnN_db/
sbatch \
  -D $PWD \
  --output $PWD/metaphlan-%j.out \
  --export=ALL \
  -J metaphlan_db_build \
  -c 8 \
  -p your_partition \
  --account=your_account \
  -t 300 \
  --wrap="metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db $INSTALL_LOCATION/metaphlan/"
```
This ensures compatibility with HUMAnN3 and prevents automatic downloads of incompatible databases.

## 5 - Phylophlan

PhyloPhlAn enables phylogenetic analysis of microbial genomes. The database setup involves downloading metadata and extracting the relevant files.

To begin, create a directory and download the necessary files:

```bash 
  mkdir phylophlan_db
  cd phylophlan_db
  wget http://cmprod1.cibio.unitn.it/databases/PhyloPhlAn/phylophlan_databases.txt
  wget $(grep -o 'http://[^[:space:]]*\.tar\b' phylophlan_databases.txt)
  wget $(grep -o 'http://[^[:space:]]*\.md5\b' phylophlan_databases.txt)
```
Next, unpack the database archives and decompress the required data files. For this pipeline, only the phylophlan db is used:
```bash
  tar -xvpf phylophlan.tar
  cd pylophlan
  bzip2 -dk phylophlan.faa.bz2
```
Then, index the databases using diamond with the following commands:
```bash
  conda activate diamond
  diamond makedb --in phylophlan/phylophlan.faa --db phylophlan/phylophlan
```
Finally, rearrange files to get a directory structure like this:
```bash
  phylophlan_db
  ├── phylophlan
  │   ├── phylophlan.dmnd
  │   ├── phylophlan.faa
  │   └── phylophlan.faa.bz2
```
In your `.env` file, `PHYLO_DB` must point to `/absolute/path/to/phylophlan_db` for PhyloPhlAn to function properly.


## 6 - GTDB

[GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) is a toolkit for assigning objective taxonomic classifications to bacterial and archaeal genomes based on the Genome Taxonomy Database (GTDB). Before using GTDB-Tk, you need to download the reference data.

Follow these steps to set up the GTDB database:

```bash 
# Manually download the latest reference data:
wget https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz

# Extract the archive to a target directory:
tar -xvzf gtdbtk_r226_data.tar.gz -C "/path/to/gtdbtk_db" --strip 1 > /dev/null
rm gtdbtk_r226_data.tar.gz
```
Verify your installation:

```bash
conda activate gtdbtk
#Set the environment variable to make GTDB-Tk aware of the database location:
export GTDBTK_DATA_PATH=/path/to/gtdbtk_db 
#verify:
gtdbtk --data-dir $GTDBTK_DATA_PATH identify 
```

## 7 - Kraken2

Kraken2 is a taxonomic classification system that assigns taxonomic labels to short DNA sequences. It uses exact alignment of k-mers and a compact hash-based index.

Kraken2 provides a [standard prebuilt database (60 GB)](https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20250402.tar.gz) for general use and many others available [here](https://benlangmead.github.io/aws-indexes/k2). However, for best results, the gtdb database formatted for kraken2 is recommended. The procedure to prepare a Kraken2 database with GTDB data is documented [here](./kraken2_gtdb.md).

