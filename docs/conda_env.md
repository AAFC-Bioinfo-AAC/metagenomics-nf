# Setting up conda environments for databases preparation

Conda environments can be used instead of apptainer containers to build the necessary databases. Here we describe the approaches to create these conda environments.  

Before starting, be sure to add the following channels to your conda installation:  

```shell
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels biobakery
```

### Biobakery3 (Chocoplhan, uniref and metaphlan databases)

```shell
conda create --name biobakery3 -c conda-forge -c bioconda -c biobakery python=3.7 humann=3.6 metaphlan=4.0.3
```

### CheckM2

```shell
conda create -n checkm2 -c bioconda checkm2=1.1.0 python=3.12
```

### DRAM

```shell
conda create -n DRAM -c bioconda dram=1.5.0 mmseqs2=17.b804f
```

### GTDB

```shell
conda create -n gtdbtk -c conda-forge -c bioconda gtdbtk=2.4.1
```

### Kaiju

```shell
conda create -n kaiju -c bioconda kaiju
```

### Kraken2/Bracken

```shell
conda create -n kraken2 -c bioconda kraken2=2.14
```

### diamond (required for phylophlan)

```shell
conda create -n diamond -c bioconda diamond
```
