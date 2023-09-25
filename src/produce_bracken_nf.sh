#!/bin/bash


# Usage: ./produce_braken.sh (i) kraken2_path (ii) baseDir as in Nextflow context
# Source: https://www.nicholas-ollberding.com/post/taxonomic-and-functional-profiling-using-biobakery-workflows
# Author : Xavier Monger, adapted by Jean-Simon Brouard
# Date : 31 aout 2021


#Variables

KRAKEN_DB=$1
# baseDir is a nextflow variable that we want to use 
baseDir=$2
echo "$baseDir:" > tester
echo "$KRAKEN_DB" >> tester
# prerequisite: avoir créé le répertoire k2_assembly_reports contenant tous les fichiers kreport
# (créés en ajoutant --report lors de l'appel à kraken)

cd k2_assembly_reports
mkdir braken
mkdir braken/species
mkdir braken/genus
mkdir braken/phylum

#dépendemment du noms données aux fichiers output de kraken, il va peut-être falloir changer les noms si *_report.txt n'est pas spécifique au fichiers désirés
#Cette étape crée des rapport bracken à partir de fichier kreport, bracken doit préalablement être installé
#Les étapes sont répétés pour les rang taxonmique désiré, un rapport est produit pour chaque rang taxonomique

for i in *report.txt
do
  filename=$(basename "$i")
  fname="${filename%report.txt}"
  bracken -d $KRAKEN_DB -i $i -r 150 -l S -o ${fname}report_species.txt
done
mv *_bracken*.txt braken/species/.

for i in *report.txt
do
  filename=$(basename "$i")
  fname="${filename%report.txt}"
  bracken -d $KRAKEN_DB -i $i -r 150 -l G -o ${fname}report_genus.txt
done
rm *_genus.txt
mv *_bracken*.txt braken/genus/.

for i in *report.txt
do
  filename=$(basename "$i")
  fname="${filename%report.txt}"
  bracken -d $KRAKEN_DB -i $i -r 150 -l P -o ${fname}report_phylum.txt
done
rm *_phylum.txt
mv *_bracken*.txt braken/phylum/.

#Ici, les fichiers bracken sont convertis en fichier MPA et puis les rapport pour tous les échantillons sont combinés, produisant ainsi un tableau d'abondance regroupant tous les échantillons

mkdir braken/species/mpa
mkdir braken/genus/mpa
mkdir braken/phylum/mpa

for i in braken/species/*bracken*.txt
do
  filename=$(basename "$i")
  fname="${filename%report_bracken.txt}"
  $baseDir/kreport2mpa.py -r $i -o braken/species/mpa/${fname}_mpa.txt --display-header
done

mkdir braken/species/mpa/combined
$baseDir/combine_mpa.py -i braken/species/mpa/*_mpa.txt -o braken/species/mpa/combined/combined_species_mpa.txt
grep -E "(s__)|(#Classification)" braken/species/mpa/combined/combined_species_mpa.txt > braken/species/mpa/combined/bracken_abundance_species_mpa.txt


for i in braken/genus/*bracken*.txt
do
  filename=$(basename "$i")
  fname="${filename%report_bracken.txt}"
  $baseDir/kreport2mpa.py -r $i -o braken/genus/mpa/${fname}_mpa.txt --display-header
done

mkdir braken/genus/mpa/combined
$baseDir/combine_mpa.py -i braken/genus/mpa/*_mpa.txt -o braken/genus/mpa/combined/combined_genus_mpa.txt
grep -E "(g__)|(#Classification)" braken/genus/mpa/combined/combined_genus_mpa.txt > braken/genus/mpa/combined/bracken_abundance_genus_mpa.txt


for i in braken/phylum/*bracken*.txt
do
  filename=$(basename "$i")
  fname="${filename%report_bracken.txt}"
  $baseDir/kreport2mpa.py -r $i -o braken/phylum/mpa/${fname}_mpa.txt --display-header
done

mkdir braken/phylum/mpa/combined
$baseDir/combine_mpa.py -i braken/phylum/mpa/*_mpa.txt -o braken/phylum/mpa/combined/combined_phylum_mpa.txt
grep -E "(p__)|(#Classification)" braken/phylum/mpa/combined/combined_phylum_mpa.txt > braken/phylum/mpa/combined/bracken_abundance_phylum_mpa.txt

mkdir ../bracken_abundance_files

#Cleaning up sample names

sed  -r -e 's/Kraken2_(\S+)\.report_bracken_[a-z]+.txt/\1/g' braken/species/mpa/combined/bracken_abundance_species_mpa.txt > ../bracken_abundance_files/bracken_abundance_species_mpa.txt
sed  -r -e 's/Kraken2_(\S+)\.report_bracken_[a-z]+.txt/\1/g' braken/genus/mpa/combined/bracken_abundance_genus_mpa.txt > ../bracken_abundance_files/bracken_abundance_genus_mpa.txt
sed  -r -e 's/Kraken2_(\S+)\.report_bracken_[a-z]+.txt/\1/g' braken/phylum/mpa/combined/bracken_abundance_phylum_mpa.txt > ../bracken_abundance_files/bracken_abundance_phylum_mpa.txt

cp -r braken ..

# original
#sed -i -e 's/report_bracken.txt//g' braken/species/mpa/combined/bracken_abundance_species_mpa.txt
#sed -i -e 's/report_bracken.txt//g' braken/genus/mpa/combined/bracken_abundance_genus_mpa.txt
#sed -i -e 's/report_bracken.txt//g' braken/phylum/mpa/combined/bracken_abundance_phylum_mpa.txt





