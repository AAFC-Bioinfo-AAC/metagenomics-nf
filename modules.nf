process QUALITY_FILTERING {

  publishDir "$projectDir/fastp"

  //label "mem_small" 
  input: 
    tuple val(datasetID), path(read1), path(read2)
 
  output:   
    tuple \
      val(datasetID), \
      path("${datasetID}_trimmed_R1.fastq.gz"), \
      path("${datasetID}_trimmed_R2.fastq.gz"), \
      path("${datasetID}_unpaired_R1.fastq.gz"), \
      path("${datasetID}_unpaired_R2.fastq.gz"), \
      path("${datasetID}.html")
  
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



process BOWTIE2 {

  publishDir "$projectDir/bowtie2"

  //label "mem_large" 
  input:
    path genome, stageAs: "genome" //to rename the folder containing the bt2 files
    val(genome_basename)
    tuple \
      val(datasetID), \
      path(trimmed_R1), \
      path(trimmed_R2), \
      path(unpaired_R1), \
      path(unpaired_R2), \
      path(html)
 
  output:   
    tuple \
      val(datasetID), \
      path("${datasetID}.sam")
  
  script:
  """
  bowtie2 -x "genome/${genome_basename}" -1 ${trimmed_R1} -2 ${trimmed_R2} -S ${datasetID}.sam
  """
}


process SAM2BAM {

  publishDir "$projectDir/unmapped_bam"

  input: 
    tuple \
      val(datasetID), \
      path(aln)
 
  output:   
    tuple \
      val(datasetID), \
      path("${datasetID}_unmapped.bam")
  
  script:
  """
  samtools view -@ 30 -bS -f 12 -F 256 ${aln} > ${datasetID}_unmapped.bam
  """
}

process SORTBAM {

  publishDir "$projectDir/unmapped_sorted_bam"

  input: 
    tuple \
      val(datasetID), \
      path(aln)
 
  output:   
    tuple \
      val(datasetID), \
      path("${datasetID}_unmapped.sorted.bam")
  
  script:
  """
  samtools sort -@ 2 -n ${aln} > ${datasetID}_unmapped.sorted.bam
  """
}


process OUTPUT_UNALIGNED_READS {
  publishDir "$projectDir/unaligned_reads"

  input: 
    tuple \
      val(datasetID), \
      path(aln)
 
  output:   
    tuple \
      val(datasetID), \
      path("${datasetID}_R1.fastq.gz"), \
      path("${datasetID}_R2.fastq.gz")
  
  script:
  """
  bedtools bamtofastq -i ${aln} -fq ${datasetID}_R1.fastq -fq2 ${datasetID}_R2.fastq
  gzip ${datasetID}_R1.fastq
  gzip ${datasetID}_R2.fastq
  """
}


process KAIJU {
  publishDir "$projectDir/kaiju_1"

  input:
    path db
    tuple \
      val(datasetID), \
      path(final_R1), \
      path(final_R2)
 
  output:   
    tuple \
      val(datasetID), \
      path("${datasetID}.out")
  
  script:
  """
  kaiju -t ${db}/nodes.dmp \
        -f ${db}/kaiju_db_nr_euk.fmi \
        -i ${final_R1} \
        -j ${final_R2} \
        -v \
        -o ${datasetID}.out \
        -z 40
  """
}


process KAIJU_TAX_TABLE {
  publishDir "$projectDir/kaiju_tax_table"

  input:
    path db
    tuple \
      val(datasetID), \
      path(kaiju_out)
 
  output:   
    tuple \
      val(datasetID), \
      path("${datasetID}.species.summary.tsv")
  
  script:
  """
  kaiju2table -t ${db}/nodes.dmp \
        -n ${db}/names.dmp \
        -r species \
        -o ${datasetID}.species.summary.tsv \
        ${kaiju_out}
  """
}


process KAIJU_FULL_TAX_TABLE {
  publishDir "$projectDir/kaiju_full_tax_table"

  input:
    path db
    tuple \
      val(datasetID), \
      path(kaiju_out)
 
  output:   
    tuple \
      val(datasetID), \
      path("${datasetID}.all_tax.summary.tsv")
  
  script:
  """
  kaiju2table -t ${db}/nodes.dmp \
        -n ${db}/names.dmp \
        -r species \
        -l superkingdom,phylum,class,order,family,genus,species \
        -o ${datasetID}.all_tax.summary.tsv \
        ${kaiju_out}
  """
}



process MERGE_TAX_FILES {
  publishDir "$projectDir/kaiju_merged"

  input:
      path (tsv_files, stageAs: "Species/*")   // to put input files in a folder parsed by R script

  output:   
    path("kaiju_merged_species.csv")
  
  script:
  """
  merge_tax_files.R
  """
}

process CAT_FASTQ {
  input:
    tuple \
      val(datasetID), \
      path(final_R1), \
      path(final_R2)

  output:
    tuple \
      val(datasetID), \
      path("${datasetID}_cat.fastq.gz")
    
  script:
  """
  cat ${final_R1} ${final_R2} > ${datasetID}_cat.fastq.gz
  """
}


process HUMANN_RUN {

  publishDir "$projectDir/humann_run"
  input:
    path chocophlan_db
    path metaphlan_db
    path uniref_db
    tuple \
      val(datasetID), \
      path(reads)

  output:
    tuple \
      val(datasetID), \
      path("${datasetID}_humann3_output")

  script:
  """
  humann -i $reads -o ${datasetID}_humann3_output \
         --threads 40 \
         --remove-temp-output \
         --metaphlan-options "--bowtie2db ${metaphlan_db} --index mpa_vJan21_CHOCOPhlAnSGB_202103" \
         --nucleotide-database $chocophlan_db \
         --protein-database $uniref_db
  """
}


process HUMANN_ABUNDANCE {

publishDir "$projectDir/humann_results"

input:
      path (humann3_output, stageAs: "pathway_abundance_files/*")

  output:
      path ("Joined_pathabundance.tsv")
      path ("Joined_pathabundance_relab.tsv")
      path ("Split_humann_pathwayabundance_tables/*")

  script:
  """
  mkdir tsv
  find -L . -iname "*.tsv" | xargs -i cp {} ./tsv # cp tsv files nested in sub-folders..
  humann_join_tables -i tsv -o Joined_pathabundance.tsv --file_name pathabundance
  
  # Then normalize to relative abundance values
  humann_renorm_table \
        --input Joined_pathabundance.tsv \
        --output Joined_pathabundance_relab.tsv \
        --units relab 

  # Finally split the Joinded_pathabundance.tsv table into "stratified" and
  # "unstratified" tables.

  mkdir Split_humann_pathwayabundance_tables

  humann_split_stratified_table \
    --input Joined_pathabundance_relab.tsv \
    --output Split_humann_pathwayabundance_tables
  """
}


process COASSEMBLY {

label 'megahit'
publishDir "$projectDir/coassembly_input"

input:
  path (readsR1, stageAs: "readsR1/*")
  path (readsR2, stageAs: "readsR2/*")

output:

  path ("coassembly_R1.fastq.gz")
  path ("coassembly_R1.fastq.gz")
  path ("Megahit_coassembly/*")


script:
"""
cat $readsR1 > coassembly_R1.fastq.gz
cat $readsR2 > coassembly_R2.fastq.gz

megahit -1 coassembly_R1.fastq.gz -2 coassembly_R2.fastq.gz -o Megahit_coassembly --out-prefix Coassembly -t 30 --min-contig-len 1000
"""

}

