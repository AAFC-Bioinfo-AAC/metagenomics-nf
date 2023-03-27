# metagenomic_nf


## 1 - Introduction

Write an introduction here...



## 2 - Install Nextflow

This pipeline is build with the Nextflow language. If you are not familiar with this, you can read more here.



## 3 - Set up conda environments

Nextflow and conda environments work well together on the biocluster.

Almost all process of the workflow use specific conda environment. For the sake of code portability, our approach is to declare in the nextflow.config file the path of the conda environment required for each process.

## 3.1 - Pre-build conda environments exist on the biocluster

On the biocluster, you can use my conda environemnts or copy them :

```shell
my path
```


## 3.2 - Pre-build conda environments exist on the biocluster

Here we described the command lines that have been used to create some of the conda environments.

### R

```shell
conda create -n R -c conda-forge r-base=4.2.3 r-tidyr=1.3.0
```

### biobakery3 (Humann3)

According to instructions given [here](https://huttenhower.sph.harvard.edu/humann).

First be sure to add the following channels :

```shell
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels biobakery
```

The following conda command works for me:

```shell
conda create --name biobakery3 -c conda-forge -c bioconda -c biobakery python=3.7 humann=3.6 metaphlan=4.0.3
```


## Biocluster launch :rocket:

```shell
 nextflow run main.nf -c nextflow.config -profile biocluster -resume --with report my_report
```

## Local launch (qcshera684498 miniserver)

```shell
 nextflow run main.nf -c nextflow.config -profile local -resume --with report my_report
```




## 4 - Humann3 preparation

### 4.1 - Download databases with the built in humann_databases command

Once you have a working conda environment, you will need to download some databases.

Define a standard location for these databases, here we put them in the project folder :


#### 4.1.1 - Pangenome database

```shell
screen -S pangenome
conda activate biobakery3
INSTALL_LOCATION=$PWD/db
#To upgrade your pangenome database: 
humann_databases --download chocophlan full $INSTALL_LOCATION --update-config yes

ctlr + A + D (exit screen)
```

#### 4.1.2 - Protein database

```shell
screen -S uniref
conda activate biobakery3
INSTALL_LOCATION=$PWD/db
#To upgrade your protein database: 
humann_databases --download uniref uniref90_diamond $INSTALL_LOCATION --update-config yes

ctlr + A + D (exit screen)
```

#### 4.1.3 - Annotation database

```shell
screen -S utility
conda activate biobakery3
INSTALL_LOCATION=$PWD/db
#To upgrade your ut database: 
humann_databases --download utility_mapping full $INSTALL_LOCATION --update-config yes

ctlr + A + D (exit screen)
```


### 4.2 - Verify your installation

You can check the location of these databases with :

```shell
humann_config --print
```

> Your interpretation is correct. HUMAnN 3.6 includes compatibility files to associate MetaPhlAn 4 profiles from the Jan21 database with the HUMAnN 3 pangenome collection. We are working on a HUMAnN 3.7 release that will add compatibility with the new database. In the meantime you can add --metaphlan-options="--index mpa_vJan21_CHOCOPhlAnSGB_202103" to your HUMAnN 3.6 command to tell MetaPhlAn 4 to get/use the previous index.





Note that despite its name, *mpa_vJan21_CHOCOPhlAnSGB_202103* appears to be the last version!

<p align="left"><img src="misc/index_chocophlan_db.png"></p>

I've also encountered this error message :

```shell
CRITICAL ERROR: The directory provided for ChocoPhlAn contains files ( mpa_vOct22_CHOCOPhlAnSGB_202212.1.bt2l ) that are not of the expected version. Please install the latest version of the database: v201901_v31
```



### How to fix an annoying bug that lead MetaPhlAn from trying to download the latest database!


#### 1 - Understand that others encounterd this issue

> Thank you for posting the command! If you install the database at the default location or if you install it in a custom location and add --bowtiedb <path/to/metaphlan/database> to your --metaphlan-options it should stop MetaPhlAn from trying to download the latest database.

> Thank you very much for the quick and helpful reply. Downgrading Metaphlan to v 4.0.3 really got rid of the initial problem. Anyhow, I still can’t stop Metaphlan from downloading the new database when I run it, even when I safe the “v21” database in the right folder. Do you know how to stop it from doing so?


  * https://forum.biobakery.org/t/humann3-compatibility-problem-with-metaphlan4/4894/2

  * https://forum.biobakery.org/t/humann3-compatibility-problem-with-metaphlan4/4894/5


#### 2 - Install the **good** database








To install the latest compatible db i *mpa_vJan21_CHOCOPhlAnSGB_202103* at a specific location, do :

```shell
metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db /isilon/sherbrooke-rdc/users/brouardjs/metagenomic_nf/data/humann3/metaphlan_vJan21/

Downloading MetaPhlAn database
Please note due to the size this might take a few minutes

\Downloading and uncompressing indexes

Downloading http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vJan21_CHOCOPhlAnSGB_202103_bt2.tar
Downloading file of size: 19909.64 MB
19909.64 MB 100.00 %  24.52 MB/sec  0 min -0 sec
Downloading http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vJan21_CHOCOPhlAnSGB_202103_bt2.md5
Downloading file of size: 0.00 MB
0.01 MB 11070.27 %  34.68 MB/sec  0 min -0 sec
Downloading and uncompressing additional files

Downloading http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.tar
Downloading file of size: 2623.07 MB
2623.07 MB 100.00 %  24.74 MB/sec  0 min -0 sec
Downloading http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.md5
Downloading file of size: 0.00 MB
0.01 MB 11702.86 %  38.42 MB/sec  0 min -0 sec

Decompressing /isilon/sherbrooke-rdc/users/brouardjs/metagenomic_nf/data/humann3/metaphlan_vJan21/mpa_vJan21_CHOCOPhlAnSGB_202103_SGB.fna.bz2 into /isilon/sherbrooke-rdc/users/brouardjs/metagenomic_nf/data/humann3/metaphlan_vJan21/mpa_vJan21_CHOCOPhlAnSGB_202103_SGB.fna

Decompressing /isilon/sherbrooke-rdc/users/brouardjs/metagenomic_nf/data/humann3/metaphlan_vJan21/mpa_vJan21_CHOCOPhlAnSGB_202103_VSG.fna.bz2 into /isilon/sherbrooke-rdc/users/brouardjs/metagenomic_nf/data/humann3/metaphlan_vJan21/mpa_vJan21_CHOCOPhlAnSGB_202103_VSG.fna
Removing uncompressed databases

Download complete
The database is installed
```

### How to fix another annoying bug : NameError: name ‘metaphlan_v4_db_version’ is not defined

Just avoid MetaPhlan version 4.0.5 and downgrade to 4.0.3.

```shell
conda activate biobakery3
conda install --name biobakery3 -c conda-forge -c bioconda -c biobakery metaphlan=4.0.3
```




## Kaiju preparation


Add the merge_tax_files to your R conda env:

```shell
cp merge_tax_files.R ../envs/R/bin/
```













