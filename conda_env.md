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

### Abricate

```shell
conda create -n abricate -c conda-forge -c bioconda -c defaults abricate
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

Do not forget the pip install step!

```shell
mamba create --prefix /isilon/common/conda/brouardjs/checkm2 -c bioconda -c conda-forge checkm2
conda activate /isilon/common/conda/brouardjs/checkm2
pip install CheckM2
```


### Coverm
```shell
mamba create -n coverm_0.6.1 -c bioconda coverm=0.6.1
```


### DRAM

```shell
conda create -n DRAM_2023 -c bioconda dram=1.4.6 mmseqs2
conda activate DRAM_2023

# Very important to have the latest version of mmseqs2!!
conda install -c bioconda mmseqs2
```

Also you may need to be sure that all users can read the file TRNAinf.cm




** On NML Waffles **

I encontered permissions errors on waffles because some files are owned by root and belong to the wheel group. For example, in spite that the file exists, tRNASCAN-SE produced errors like this :

FATAL: Unable to find /Drives/O/USERS/jsbrouard/.conda/envs/DRAM_2023/bin/covels-SE executable

When examining these files, one can notice that they have no rwx permissions for the other group like this :

```shell
-rwxrwx---. 2 root      wheel   66552  2 aoû 01:38 TRNAinf-bact-SeC.cm
-rwxrwxr--. 1 jsbrouard wheel   62945 31 aoû 12:48 TRNAinf.cm
-rwxrwx---. 2 root      wheel   62677  2 aoû 01:38 TRNAinf-euk.cm
```

A workaround  is to copy the files in another location. After that, you become the owner of these files and you are able to change permissions.

I have modified two repositories in the DRAM env :
    - the tRNASCAN-SE folder
    - the bin folder

I have someting like this :


```shell
# 1 - For the tRNASCAN folder

cd /Drives/O/USERS/jsbrouard/.conda/envs/DRAM_2023/lib
# copy the whole folder to another place
cp tRNAscan-SE /Drives/O/USERS/jsbrouard/AAFC-AAC/tRNAscan-SE
rm -rf tRNAscan-SE
ln -s /Drives/O/USERS/jsbrouard/AAFC-AAC/tRNAscan-SE tRNAscan-SE

chmod -R a+rwx tRNAscan-SE/


# 2 - For all DRAM bins!
cd /Drives/O/USERS/jsbrouard/.conda/envs/DRAM_2023/bin
mkdir new_bin
cp * new_bin
mv new_bin ..
rm -rf *
mv ../new_bin/* .
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
mamba create -n kraken2 -c bioconda kraken2=2.1.2 bracken
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

### geNomad

https://github.com/apcamargo/genomad

# Create a conda environment for geNomad
conda create -n genomad -c conda-forge -c bioconda genomad









## Novel instructions : build on share folder (NML Waffles)



mamba create --prefix $share/phylophlan -c bioconda phylophlan=3.0.3

```shell
share='/Drives/O/GRDI-AMR2/share/conda/envs'

#phylophlan
sbatch -J phylophlan_conda -c 4 --mem 32G -p NMLResearch --wrap="mamba create --prefix $share/phylophlan_3.0.3 -c bioconda phylophlan=3.0.3 -y"

#R
sbatch -J R_conda -c 4 --mem 32G -p NMLResearch --wrap="mamba create --prefix $share/R -c conda-forge r-base=4.2.3 r-tidyr=1.3.0 -y"


# quast
sbatch -J quast_conda -c 4 --mem 32G -p NMLResearch --wrap="mamba create --prefix $share/quast-5.2.0 -c bioconda quast=5.2.0 -y"

#metabat2

sbatch -J quast_metabat2 -c 4 --mem 32G -p NMLResearch --wrap="mamba create --prefix $share/metabat2_2.15 -c bioconda metabat2=2.15 -y"

#megahit

sbatch -J megahit_conda -c 4 --mem 32G -p NMLResearch --wrap="mamba create --prefix $share/megahit_1.2 -c bioconda megahit=1.2 -y"


#kaiju
sbatch -J kaiju_conda -c 4 --mem 32G -p NMLResearch --wrap="mamba create --prefix $share/kaiju_1.9.2 -c bioconda kaiju=1.9.2 -y"


#gtdb
sbatch -J gtdb_conda -c 4 --mem 32G -p NMLResearch --wrap="mamba create --prefix $share/gtdbtk-2.1.1 -c conda-forge -c bioconda gtdbtk=2.1.1 -y"
conda activate $share/gtdbtk-2.1.1
# downgrade numpy 1.24
python -m pip uninstall numpy
# Reinstall numpy
python -m pip install numpy==1.23.1


#fastp
sbatch -J fastp_conda -c 4 --mem 32G -p NMLResearch --wrap="mamba create --prefix $share/fastp_0.23.4 -c bioconda fastp=0.23.4 -y"

#drep
sbatch -J drep_conda -c 4 --mem 32G -p NMLResearch --wrap="mamba create --prefix $share/drep_3.4.5 -c bioconda checkm-genome=1.2.2 drep=3.4.5 -y"


#DRAM
sbatch -J dram_conda -c 4 --mem 32G -p NMLResearch --wrap="mamba create --prefix $share/DRAM_2023 -c bioconda dram=1.4.6 mmseqs2 -y"


#coverM
sbatch -J coverm_conda -c 4 --mem 32G -p NMLResearch --wrap="mamba create --prefix $share/coverm_0.6.1 -c bioconda coverm=0.6.1 -y"

#checkm2
sbatch -J checkm2_conda -c 4 --mem 32G -p NMLResearch --wrap="mamba create --prefix $share/checkm2_1.0.1 -c bioconda -c conda-forge checkm2=1.0.1 -y"

#Do not forget the pip install step!
conda activate $share/checkm2_1.0.1
pip install --user CheckM2

#bowtie2
sbatch -J bowtie2_conda -c 4 --mem 32G -p NMLResearch --wrap="mamba create --prefix $share/bowtie2_2.5.2 -c bioconda bowtie2=2.5.2 samtools=1.18 -y"


#biobakery
sbatch -J biobakery3_conda -c 4 --mem 32G -p NMLResearch --wrap="mamba create --prefix $share/biobakery3 -c conda-forge -c bioconda -c biobakery python=3.7 humann=3.6 metaphlan=4.0.3 -y"

#bedtools

sbatch -J bedtools_conda -c 4 --mem 32G -p NMLResearch --wrap="mamba create --prefix $share/bedtools_2.31.1 -c bioconda bedtools=2.31.1 -y"




```


