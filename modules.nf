

process QUALITY_FILTERING {
  publishDir "$projectDir/fastp"

  label "mem_xlarge" 
  input: 
    tuple val(datasetID), path(read1), path(read2)
 
  output:   
    tuple \
      val(datasetID), \
      path("${datasetID}_trimmed_R1.fastq.gz"), \
      path("${datasetID}_trimmed_R2.fastq.gz"), \
      path("${datasetID}_unpaired_R1.fastq.gz"), \
      path("${datasetID}_unpaired_R2.fastq.gz")
  
  script:
  """
    fastp -i $read1 \
        -I $read2 \
        -o ${datasetID}_trimmed_R1.fastq.gz \
        -O ${datasetID}_trimmed_R2.fastq.gz \
        --unpaired1 ${datasetID}_unpaired_R1.fastq.gz \
        --unpaired2 ${datasetID}_unpaired_R2.fastq.gz \
        --length_required 100 \
        --cut_tail --cut_front \
        --cut_mean_quality 15 \
        --qualified_quality_phred 15 \
        --cut_window_size 4 \
        --detect_adapter_for_pe \
        --report_title="${datasetID}" \
        --html=${datasetID}.html
  """
}


/*
 * Process 1A: Create a FASTA genome index (.fai) with samtools for GATK
 */

process PREPARE_GENOME_SAMTOOLS { 
  tag "$genome.baseName"
  label "mem_xlarge" 
  input: 
    path genome
 
  output: 
    path "${genome}.fai"
  
  script:
  """
  samtools faidx ${genome}
  """
}


/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

process PREPARE_GENOME_PICARD {
  tag "$genome.baseName"
  label "mem_xlarge"

  input:
    path genome
  output:
    path "${genome.baseName}.dict"

  script:
  """
  gatk CreateSequenceDictionary -R $genome -O ${genome.baseName}.dict
  """
}



/*
 * Process 1C: Create a file containing the filtered and recoded set of variants
 */

process PREPARE_VCF_FILE {
  tag "$variantsFile.baseName"

  input: 
    path variantsFile
    path denylisted

  output:
    tuple \
      path("${variantsFile.baseName}.filtered.recode.vcf.gz"), \
      path("${variantsFile.baseName}.filtered.recode.vcf.gz.tbi")
  
  script:  
  """
  vcftools --gzvcf $variantsFile -c \
           --exclude-bed ${denylisted} \
           --recode | bgzip -c \
           > ${variantsFile.baseName}.filtered.recode.vcf.gz
  tabix ${variantsFile.baseName}.filtered.recode.vcf.gz
  """
}


/*
 * Process 2: Mark duplicates
 *             
 */

process MARK_DUPLICATES {
  tag "$sampleId"
  label "mem_small"

  publishDir "$projectDir/results/mark_duplicates"
    
  input: 
    tuple val(sampleId), path(bam), path(index)


  output:
    tuple \
      val(sampleId), \
      path("${sampleId}.MD.bam"), \
      path("${sampleId}.MD.bam.bai")
  
  script:
  """
  picard MarkDuplicates \
  I= $bam \
  O= ${sampleId}.MD.bam \
  M= "${sampleId}_marked_dup_metrics.txt"
  
  # Index BAM files
  samtools index ${sampleId}.MD.bam
  """
}


/*
 * Process 3: Base recalibrate to detect systematic errors in base quality scores, 
 *            select unique alignments and index
 *             
 */

process GERMLINE_GATK_RECALIBRATE {
  tag "$sampleId"
  label "mem_small"
  publishDir "$projectDir/results/recal"

  input: 
    path genome
    path index
    path dict
    tuple val(sampleId), path(bam), path(index)
    tuple path(variants_file), path(variants_file_index)

  output:
    tuple \
      val(sampleId), \
      path("${sampleId}.final.uniq.bam"), \
      path("${sampleId}.final.uniq.bam.bai")
  
  script: 
  """
  # Indel Realignment and Base Recalibration
  gatk BaseRecalibrator \
          -R $genome \
          -I $bam \
          --known-sites $variants_file \
          -O final.GERMLINE.grp 
  gatk ApplyBQSR \
          -R $genome -I $bam \
          --bqsr-recal-file final.GERMLINE.grp \
          -O ${sampleId}.final.uniq.bam &&
  # Index BAM files
  samtools index ${sampleId}.final.uniq.bam &&
  # Remove intermediate files
  readlink -f $bam | xargs rm
  """
}


/*
 * Process 4: Call variants Per-Sample (GVCF mode)
 */

process CALLING_VARIANTS_PER_SAMPLE {
  tag "$sampleId"
  label "mem_small"

  publishDir "$projectDir/results/GVCF"

  publishDir "$projectDir/results/gvcf"

  input:
    path genome
    path index
    path dict
    tuple val(sampleId), path(bam), path(bai)
 
  output: 
    tuple val(sampleId), path("${sampleId}.g.vcf.gz")

  script:
  """
	gatk --java-options "-Xmx4g" HaplotypeCaller \
	--dont-use-soft-clipped-bases true \
	-R $genome \
	-I $bam \
	-O "${sampleId}.g.vcf.gz" \
	-ERC GVCF &&

	# Remove intermediate files
    readlink -f $bam | xargs rm
  """
}


