#!/bin/bash

################################################################################
# Job array to remove host and phiX sequence contamination from trimmed fastq files 
#
# This script may be run from within a newly created folder called cleaned, for 
# example (to be defined within the script as OUTDIR). A directory structure will
# be created to organise the different kinds of output files generated.
################################################################################
#
# For job array template, refer:
# https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/script-array-job/
#
# 
# This script does the following: bowtie2 -> samtools -> bedtools -> pigz
# - map trimmed fastq to host+phiX reference
# - write out unmapped reads to fastq files 
# - gzip fastq files
#
# Installed all required sofware to a conda env called : 
# $ conda create --name env-remove_host_phiX \
#                 --channel bioconda --channel conda-forge \
#                 bowtie2 samtools bedtools pigz \
#                 --yes
#
# This is optional and only needed if you want to use reformat.sh to verify read
# order for testing purposes when samtools sort isnt used. This takes jsut a minute
# to run and can save ~20 min. time needed for samtools sort to run
# $ conda install --name env-remove_host_phiX \
#                 --channel bioconda bbmap \
#                 --yes
# 
# To use reformat.sh vpair to verify that reads appear to be correctly paired,
# based on their names, refer:
# https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/reformat-guide/

################################################################################



#########################################################################
## Section 1: Header                                                    #
#########################################################################

# Specify name to be used to identify this run
#$ -N remove_host_phiX

# Email address (change to yours)
#$ -M sara.ricci@agr.gc.ca

# Specify mailing options: b=beginning, e=end, s=suspended, n=never, a=aborted or suspended
#$ -m a

# This sets the task range in the array and step size 
#$ -t 1-18:1

### Total of 18 tasks calculated as follows (one task per sample pair):
### $ ls /isilon/projects/J-002460_LLQ/riccis/riccis_all_samples/paired_trimmed/NS*_R1*gz | wc -l

# This sets the maximum number of concurrent tasks 
#$ -tc 18

# Specify the number of processors for the single job
#$ -pe smp 5

# Change directory to the current 
#$ -cwd

# Specify that bash shell should be used to process this script
#$ -S /bin/bash

# Specify the output and err files
#$ -o $JOB_NAME_$TASK_ID.out
#$ -e $JOB_NAME_$TASK_ID.err

#########################################################################
## Section 2: Definitions                                               #
#########################################################################

# Define project dir
PROJECT_DIR="/isilon/projects/J-002460_LLQ"
  
# Define input dir
INDIR=${PROJECT_DIR}"/riccis/riccis_all_samples/paired_trimmed"

# Define list of trimmed fastq files
LIST_INPUT_FILES=($(ls ${INDIR}/*_R1*fastq.gz))

# Specify reference database directory 
DB_DIR="/isilon/projects/J-002460_LLQ/riccis/cow/cow"

# Define destination directory to store hostfree files
OUTDIR=${PROJECT_DIR}"/riccis/riccis_all_samples/cleaned_all"

# Create destination directory structure to organise output files
mkdir -p ${OUTDIR}/{logs,reports,clean_fastq}

#########################################################################
## Section 3: Executing the program                                     #
#########################################################################

# Note time 
echo -e "\n [START TIME]: `date` \n" 
SECONDS=0

# Activate conda
source ~/miniconda3/etc/profile.d/conda.sh
conda activate env-remove_host_phiX

# Note version of software used
# if [ ! -f ${OUTDIR}/software.version ]; then
#    echo -e "\nSAMTOOLS\n" >> ${OUTDIR}/software.version
#    samtools --version >> ${OUTDIR}/software.version 2>&1
#    echo -e "\n#################################\n" >> ${OUTDIR}/software.version
#    bedtools --version >> ${OUTDIR}/software.version 2>&1
#    echo -e "\n#################################\n" >> ${OUTDIR}/software.version
#    pigz --version >> ${OUTDIR}/software.version 2>&1
#    echo -e "\n#################################\n" >> ${OUTDIR}/software.version
#    bowtie2 --version >> ${OUTDIR}/software.version 2>&1
# fi
if [ ! -f ${OUTDIR}/software.version ]; then
   {
        echo -e "\n#SAMTOOLS\n" 
        samtools --version
        echo -e "\n#bedtools\n" 
        bedtools --version 
        echo -e "\n#pigz\n" 
        pigz --version 
        echo -e "\n#bowtie\n" 
        bowtie2 --version 
   } >> ${OUTDIR}/software.version 2>&1
fi

# Extract sample name from input file name pulled according to job id
SAMPLENAME=$(basename -s _trimmed_R1.fastq.gz ${LIST_INPUT_FILES[$SGE_TASK_ID - 1]})

# Print and run command
echo -e "\n Processing sample $SAMPLENAME [array job: $SGE_TASK_ID] \n"

# Define command
cmd="bowtie2 -x $DB_DIR \
      -1 ${INDIR}/${SAMPLENAME}_trimmed_R1.fastq.gz \
      -2 ${INDIR}/${SAMPLENAME}_trimmed_R2.fastq.gz \
      --threads $NSLOTS | \

      samtools view -u -f 12 -F 256 --threads $NSLOTS | \

      samtools sort -n -m 4G --threads $NSLOTS | \

      bedtools bamtofastq -i - \
      -fq ${OUTDIR}/clean_fastq/${SAMPLENAME}_trimmed_clean_R1.fastq \
      -fq2 ${OUTDIR}/clean_fastq/${SAMPLENAME}_trimmed_clean_R2.fastq"
      
echo $cmd
eval $cmd

echo -e "\n Mapping and cleaning took: `expr ${SECONDS} / 60` min. \n"

# # Validate read order. It seems samtools sort is not needed and this step verifies
# # that the read pairs are ordered by name regardless.
# cmd="reformat.sh in=${OUTDIR}/clean_fastq/${SAMPLENAME}_trimmed_clean_R#.fastq vpair"
# echo -e "\n Validating paired end order of reads\n"
# echo $cmd
# eval $cmd


# gzip fastq
cmd="pigz --quiet --processes $NSLOTS ${OUTDIR}/clean_fastq/${SAMPLENAME}*"
echo $cmd
eval $cmd

echo -e "\n [END TIME]: `date`"

echo -e "\n [RUN TIME]: `expr ${SECONDS} / 60` min."

# Move log files to folder called logs
sleep 2
mv *_$SGE_TASK_ID.out ${OUTDIR}/logs/${SAMPLENAME}.log
mv *_$SGE_TASK_ID.err ${OUTDIR}/reports/${SAMPLENAME}_mapping.stats


