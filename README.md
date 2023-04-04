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
Here we describe the command lines that have been used to create some of the conda environments.

### R
```shell
conda create -n R -c conda-forge r-base=4.2.3 r-tidyr=1.3.0
```


### biobakery3 (Humann3)
The instructions are given [here](https://huttenhower.sph.harvard.edu/humann).

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

### megahit
```shell
conda create -n megahit -c bioconda megahit=1.2.
```


## 4 - Humann3 preparation
### 4.1 - Download databases with the built in humann_databases command
Once you have a working conda environment, you will need to download some databases.

Define a standard location for these databases, here we put them in the project folder :

#### 4.1.1 - Pangenome database
```shell
screen -S pangenome
conda activate biobakery3
INSTALL_LOCATION=$PWD/db/humann3
#To upgrade your pangenome database: 
humann_databases --download chocophlan full $INSTALL_LOCATION --update-config yes

ctlr + A + D (exit screen)
```

#### 4.1.2 - Protein database
```shell
screen -S uniref
conda activate biobakery3
INSTALL_LOCATION=$PWD/db/humann3
#To upgrade your protein database: 
humann_databases --download uniref uniref90_diamond $INSTALL_LOCATION --update-config yes

ctlr + A + D (exit screen)
```

#### 4.1.3 - Annotation database
```shell
screen -S utility
conda activate biobakery3
INSTALL_LOCATION=$PWD/db/humann3
#To upgrade your ut database: 
humann_databases --download utility_mapping full $INSTALL_LOCATION --update-config yes

ctlr + A + D (exit screen)
```


### 4.2 - Verify your installation
You can check the location of these databases with :

```shell
humann_config --print
```


### 4.3 - Fix an annoying bug that lead MetaPhlAn from trying to download the latest database!
As is, with Humann3.6 and Metephlan 4.0.3, you will get an error message like this  :

```shell
  Downloading MetaPhlAn database
  Please note due to the size this might take a few minutes

  File /home/brouardjs/miniconda3/envs/biobakery3/lib/python3.7/site-packages/metaphlan/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212.tar already present!

  Downloading http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212.md5
  Downloading file of size: 0.00 MB
  0.01 MB 11702.86 %  74.14 MB/sec  0 min -0 sec
  MD5 checksums do not correspond! If this happens again, you should remove the database files and rerun MetaPhlAn so they are re-downloaded
```


We investigate this issue and found that this occur because metphlan will try to download and use the very latest version of database that is not compatible with humann3 at this time!!


> Your interpretation is correct. HUMAnN 3.6 includes compatibility files to associate MetaPhlAn 4 profiles from the Jan21 database with the HUMAnN 3 pangenome collection. We are working on a HUMAnN 3.7 release that will add compatibility with the new database. In the meantime you can add --metaphlan-options="--index mpa_vJan21_CHOCOPhlAnSGB_202103" to your HUMAnN 3.6 command to tell MetaPhlAn 4 to get/use the previous index.
...

> Thank you for posting the command! If you install the database at the default location or if you install it in a custom location and add --bowtiedb <path/to/metaphlan/database> to your --metaphlan-options it should stop MetaPhlAn from trying to download the latest database.
...


> Thank you very much for the quick and helpful reply. Downgrading Metaphlan to v 4.0.3 really got rid of the initial problem. Anyhow, I still can’t stop Metaphlan from downloading the new database when I run it, even when I safe the “v21” database in the right folder. Do you know how to stop it from doing so?


  * https://forum.biobakery.org/t/humann3-compatibility-problem-with-metaphlan4/4894/2

  * https://forum.biobakery.org/t/humann3-compatibility-problem-with-metaphlan4/4894/5


The take home message is that the latest db that can be used by Metaphlan 3.6 is *mpa_vJan21_CHOCOPhlAnSGB_202103* (not mpa_vOct22_CHOCOPhlAnSGB_202212).


### 4.3.1 - Solution: install the **latest compatible db** with Metaphlan3.6
To install the latest compatible db i *mpa_vJan21_CHOCOPhlAnSGB_202103* at a specific location, do :

```shell
screen -S metaphlan
conda activate biobakery3
metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db /mnt/sdb/metagenomic_nf/db/humann3/metaphlan_vJan21/
```

And you are all set; we have added the right command lines options to specify which database should be used by Metaphlan when Humann3 is invoked!



## 5 - Kaiju preparation
Add the merge_tax_files to your R conda env:

```shell
cp -v src/merge_tax_files.R /home/brouardjs/miniconda3/envs/R/bin/
```



## 6 - Launch :rocket:

### 6.1 - Biocluster
```shell
 nextflow run main.nf -c nextflow.config -profile biocluster -resume --with report my_report
```

### 6.2 - Local server (qcshera684498 miniserver)
```shell
 nextflow run main.nf -c nextflow.config -profile local -resume --with report my_report
```











