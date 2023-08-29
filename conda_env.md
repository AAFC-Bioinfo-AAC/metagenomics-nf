# 3 - Set up conda environments

For the sake of code portability, we declare the conda environment required for each process in the nextflow.config file.

Here we describe the approaches that have been used to create these conda environments.


Before starting, be sure to add the following channels to your conda installation :

```shell
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels biobakery
```


### Bedtools
```shell
 conda create -n bedtools -c bioconda bedtools
```


### Biobakery3 (Humann3)

The instructions are given [here](https://huttenhower.sph.harvard.edu/humann).

The following conda command works for me:

```shell
conda create --name biobakery3 -c conda-forge -c bioconda -c biobakery python=3.7 humann=3.6 metaphlan=4.0.3
```

### Bowtie2

Note that the Samtools are also required in thie environment.

```shell
mamba create -n bowtie2 -c bioconda bowtie2 samtools
```

### CheckM2

This works for me :

```shell
git clone --recursive https://github.com/chklovski/checkm2.git && cd checkm2
mamba create -n checkm2 -f checkm2
conda activate checkm2
pip install CheckM2
```


### Coverm
```shell
mamba create -n coverm -c bioconda coverm
```


### DRAM

```shell
conda create -n DRAM_2023 -c bioconda dram
conda activate DRAM_2023

# Very important to have the latest version of mmseqs2!!
conda install -c bioconda mmseqs2
```

Note that the wiki has a plenty of information.

[Reference](https://github.com/WrightonLabCSU/DRAM)



### Drep
```shell
mamba create -n drep -c bioconda checkm-genome drep
```


### Fastp

```shell
 conda create -n fastp -c bioconda fastp
```


### GTDB

https://ecogenomics.github.io/GTDBTk/installing/index.html

```shell
mamba create -n gtdbtk-2.1.1 -c conda-forge -c bioconda gtdbtk=2.1.1

conda activate gtdbtk-2.1.1

# downgrade numpy 1.24
python -m pip uninstall numpy
# Reinstall numpy
python -m pip install numpy==1.23.1
```


### Kaiju

```shell
mamba create -n kaiju -c bioconda kaiju
```


### Kraken2/Bracken

Intall Bracken alongside kraken2

```shell
conda create -n kraken2 -c bioconda kraken2=2.1.2
conda activate kraken2
conda install -c bioconda bracken
```



### Megahit

```shell
conda create -n megahit -c bioconda megahit=1.2.
```


### Metabat2

```shell
conda create -n metabat2 -c bioconda metabat2=2.15
```


### Phylophlan
```shell
mamba create -n phylophlan -c bioconda phylophlan=3.0.3
```


### Quast

```shell
mamba create -n quast-5.2.0 -c bioconda quast=5.2.0
```

### R
```shell
conda create -n R -c conda-forge r-base=4.2.3 r-tidyr=1.3.0
```












