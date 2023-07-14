# Metagenomics databases preparation




## 1 - Kaiju

...

## 2 - Checkm2

```shell
conda activate checkm2
#checkm2 database --download --path /custom/path/
checkm2 database --download --path /Drives/K/jsbrouard/AAFC-AAC/db/checkm2
```

>The database path can also be set by setting the environmental variable CHECKM2DB using:

```shell
export CHECKM2DB="path/to/database"
```


## 3 - Humann3

### 3.1 - Download databases with the built in humann_databases command
Once you have a working conda environment, you will need to download some databases.

Define a standard location for these databases, here we put them in the project folder :

#### 3.1.1 - Pangenome database

```shell
screen -S pangenome
conda activate biobakery3
INSTALL_LOCATION=HUMAnN_db/HUMAnN_db_20230710
#To upgrade your pangenome database: 
humann_databases --download chocophlan full $INSTALL_LOCATION --update-config yes

ctlr + A + D (exit screen)
```

#### 3.1.2 - Protein database

```shell
screen -S uniref
conda activate biobakery3
INSTALL_LOCATION=HUMAnN_db/HUMAnN_db_20230710
#To upgrade your protein database: 
humann_databases --download uniref uniref90_diamond $INSTALL_LOCATION --update-config yes

ctlr + A + D (exit screen)
```

#### 3.1.3 - Annotation database

```shell
screen -S utility
conda activate biobakery3
INSTALL_LOCATION=HUMAnN_db/HUMAnN_db_20230710
#To upgrade your ut database: 
humann_databases --download utility_mapping full $INSTALL_LOCATION --update-config yes

ctlr + A + D (exit screen)
```


### 3.1.4 - Verify your installation
You can check the location of these databases with :

```shell
humann_config --print
```

## 4 - Metaphlan

### 4.1 - Fix an annoying bug that lead MetaPhlAn from trying to download the latest database!
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


### 4.2 - Solution: install the **latest compatible db** with Metaphlan3.6
To install the latest compatible db i *mpa_vJan21_CHOCOPhlAnSGB_202103* at a specific location, do :

```shell
screen -S metaphlan
conda activate biobakery3
metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db /mnt/sdb/metagenomic_nf/db/humann3/metaphlan_vJan21/
```

On waffles !
```shell
source activate biobakery3

sbatch -D $PWD --output $PWD/metaphlan_db_build-$j.out --export=ALL -J metaphlan_db_build -c 8 -p NMLResearch --mem 256G --wrap="metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 --bowtie2db /Drives/K/jsbrouard/AAFC-AAC/db/HUMAnN_db/HUMAnN_db_20230710/metaphlan_vJan21/"
```


And you are all set; we have added the right command lines options to specify which database should be used by Metaphlan when Humann3 is invoked!


