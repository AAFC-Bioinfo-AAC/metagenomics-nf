process RENAME_SEQUENCES {

  publishDir "$baseDir/data"

  input:
    path (seq, stageAs: "data/*")
    path map_file
    
  output:
    tuple \
      path("log.txt"), \
      path("renamed/*")
  
  script:
  """
  mkdir renamed
  $baseDir/src/rename_sequences.py $map_file \$PWD/data renamed/ fastq.gz > log.txt
  """
}

process QUALITY_FILTERING {
  label 'cpus_medium'

  publishDir "$projectDir/results/trimmed_reads"

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
        --thread=$task.cpus \
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

  label 'mem_medium'
  label 'cpus_large'
  
  publishDir "$projectDir/results/decontamination/genome/bowtie2"
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
      path("${datasetID}_unmapped.sorted.bam")
  
  script:
  """
  bowtie2 -p $task.cpus -x "genome/${genome_basename}" -1 ${trimmed_R1} -2 ${trimmed_R2} | \
  
  samtools view -u -f 12 -F 256 --threads $task.cpus | \

  samtools sort -n -m 4G --threads $task.cpus | \

  # step3 (old SORTBAM process)
  samtools sort -n -m 4G --threads $task.cpus > ${datasetID}_unmapped.sorted.bam
  """
}


process OUTPUT_UNALIGNED_READS {

  publishDir "$projectDir/results/prepared_reads"

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
  bedtools bamtofastq -i ${aln} -fq ${datasetID}_R1.fastq -fq2 ${datasetID}_R2.fastq &&
  gzip ${datasetID}_R1.fastq &&
  gzip ${datasetID}_R2.fastq
  """
}


process KAIJU {

  label 'mem_medium'
  label 'cpus_large'
  
  publishDir "$projectDir/results/kaiju/kaiju_1"

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
        -z $task.cpus
  """
}


process KAIJU_TAX_TABLE {

  label 'mem_medium'
  
  publishDir "$projectDir/results/kaiju/kaiju_tax_table"

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

  label 'mem_medium'
  
  publishDir "$projectDir/results/kaiju/kaiju_full_tax_table"

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

  publishDir "$projectDir/results/kaiju/kaiju_merged"

  input:
    path (tsv_files, stageAs: "Species/*")   // to put input files in a folder parsed by R script

  output:
    path("kaiju_merged_species.csv")
  
  script:
  """
  $baseDir/src/merge_tax_files.R
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

  label 'mem_medium'
  label 'cpus_xxlarge'

  publishDir "$projectDir/results/humann/humann_run"
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
  # Sometimes humann3 is unable to write in the output
  # e.g. :
  # 40K drwxrwx---. 2 jsbrouard grp_jsbrouard   22 25 aoû 13:10 G2-E1-13_humann3_output
  # Try to create the output folder before humann3 and change permissions
  mkdir -p ${datasetID}_humann3_output &&
  chmod a+rwx ${datasetID}_humann3_output &&

  humann -i $reads -o ${datasetID}_humann3_output \
         --threads $task.cpus \
         --remove-temp-output \
         --metaphlan-options "--bowtie2db ${metaphlan_db} --index mpa_vJan21_CHOCOPhlAnSGB_202103" \
         --nucleotide-database $chocophlan_db \
         --protein-database $uniref_db
  """
}


process HUMANN_ABUNDANCE {

  publishDir "$projectDir/results/humann/humann_results"

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


// modules related to co-assemblies

process COASSEMBLY {

  label 'mem_xxlarge'
  label 'cpus_xxlarge'
  
  publishDir "$projectDir/results/coassembly/megahit"

  input:
    path (readsR1, stageAs: "readsR1/*")
    path (readsR2, stageAs: "readsR2/*")

  output:
    path ("Megahit_coassembly/*")

  script: 
  """
  cat $readsR1 > coassembly_R1.fastq.gz
  cat $readsR2 > coassembly_R2.fastq.gz
  megahit -1 coassembly_R1.fastq.gz \
          -2 coassembly_R2.fastq.gz \
          -o Megahit_coassembly \
          --out-prefix Coassembly \
          -t $task.cpus \
          --min-contig-len 1000 \
          --presets 'meta-large'
  """
}


process BOWTIE2_BUILD {
  

  publishDir "$projectDir/results/coassembly/bwt2_index"

  input:
    path (megahit_coassembly_outfiles, stageAs: "megahit/*")

  output:
    path ("coassembly/*")
  
  script:
  """
  mkdir coassembly
  bowtie2-build megahit/Coassembly.contigs.fa coassembly/coassembly
  """
}



process BOWTIE2_MAP {
  label 'cpus_large'
  
  publishDir "$projectDir/results/coassembly/bwt2_output_for_metabat"

input:
  path bwt2_index, stageAs: "coassembly/*" //to rename the folder containing the bt2 files
  tuple \
    val(datasetID), \
    path(final_R1), \
    path(final_R2)

output:
  tuple \
    val(datasetID), \
    path("${datasetID}_sorted.bam")

script:
"""
bowtie2 -x coassembly/coassembly \
        -1 ${final_R1} \
        -2 ${final_R2} \
        -p $task.cpus | \
        samtools sort \
        --threads $task.cpus \
        -O BAM -o ${datasetID}_sorted.bam

"""
}




process JGI_SUMMARIZE {
  
  publishDir "$projectDir/results/coassembly/jgi"

  input: 
    tuple \
      val(datasetID), \
      path(aln)
 
  output:
    tuple \
      val(datasetID), \
      path("${datasetID}_depth.txt")
  
  script:
  """
  jgi_summarize_bam_contig_depths \
    --outputDepth ${datasetID}_depth.txt \
    $aln
  """
}


process METABAT2_BIN_COASSEMBLY {

  label 'cpus_xlarge'
  
  publishDir "$projectDir/results/coassembly/metabat2_bins"
  
  input:
    path (megahit_coassembly_outfiles, stageAs: "megahit/*")
    tuple \
      val(datasetID), \
      path("${datasetID}_depth.txt")
 
  output:
    tuple \
      val(datasetID), \
      path("${datasetID}/*"), optional: true
  
  script:
  """
  mkdir ${datasetID}
  metabat2 \
    -i megahit/Coassembly.contigs.fa \
    -a ${datasetID}_depth.txt \
    -o ${datasetID}/${datasetID}.bin \
    -t $task.cpus \
    -m 2000 \
    -v
  """
}


process CHECKM {

  label 'cpus_medium'
  
  publishDir "$projectDir/results/coassembly/checkM2_output"
  
  input:
    path db
    tuple \
      val(datasetID), \
      path (metabat2_coassembly_outfiles, stageAs: "Coassembled_bins/*")

  output:
    tuple \
      val(datasetID), \
      path("$datasetID/*")
  
  script:
  """
  export HDF5_USE_FILE_LOCKING='FALSE'
  checkm2 predict \
    --database_path $db \
    --threads 20 \
    -x fa \
    --input Coassembled_bins \
    --output_directory ${datasetID}
  """
}

// modules based on individual assemblies

process MEGAHIT_SINGLE {
  
  label 'mem_xxlarge'
  label 'cpus_xlarge'
  
  publishDir "$projectDir/results/indiv_assemblies/megahit"

  input:
    tuple \
      val(datasetID), \
      path(final_R1), \
      path(final_R2)

  output:
    tuple \
      val(datasetID), \
      path ("${datasetID}/*")

  script:
  """
  megahit -1 ${final_R1} \
          -2 ${final_R2} \
          -o ${datasetID} \
          --out-prefix ${datasetID} \
          -t $task.cpus \
          --min-contig-len 1000 &&

  # Remove intermediate files (intermediate contigs)
  rm -rf ${datasetID}/intermediate_contigs
  """
}



process BOWTIE2_BUILD_SINGLE {

  
  label 'cpus_xlarge'
  publishDir "$projectDir/results/indiv_assemblies/bwt2_index"

  input:
    tuple \
      val(datasetID), \
      path(megahit_individual_outfiles, stageAs: "megahit/*")

  output:
    tuple \
      val(datasetID), \
      path ("bwt2_index/*")
  
  script:
  """
  mkdir bwt2_index
  bowtie2-build \
    megahit/${datasetID}.contigs.fa \
    bwt2_index/${datasetID} \
    --quiet \
    --threads $task.cpus
  """
}


process BOWTIE2_MAP_SINGLE {

  
  label 'cpus_large'

  publishDir "$projectDir/results/indiv_assemblies/bwt2_output_for_metabat"

  input:
    tuple \
      val(datasetID), \
      path ("bwt2_index/*", stageAs: "index/*"), \
      path(final_R1), \
      path(final_R2)

  output:
    tuple \
      val(datasetID), \
      path("${datasetID}_sorted.bam")

  script:
  """
  bowtie2 -x index/${datasetID} \
          -1 ${final_R1} \
          -2 ${final_R2} \
          -p $task.cpus | \
          samtools sort \
          --threads $task.cpus \
          -O BAM -o ${datasetID}_sorted.bam
          
  """
}


process JGI_SUMMARIZE_SINGLE {

  publishDir "$projectDir/results/indiv_assemblies/jgi"

  input: 
    tuple \
      val(datasetID), \
      path(aln)

  output:
    tuple \
      val(datasetID), \
      path("${datasetID}_depth.txt")
  
  script:
  """
  jgi_summarize_bam_contig_depths \
    --outputDepth ${datasetID}_depth.txt \
    $aln
  """
}



process METABAT2_BIN_SINGLE {

  label 'cpus_xlarge'

  publishDir "$projectDir/results/indiv_assemblies/metabat2_bins"
  
  input:
    tuple \
      val(datasetID), \
      path("${datasetID}_depth.txt"), \
      path(megahit_individual_outfiles, stageAs: "megahit/*")

  output:
    tuple \
      val(datasetID), \
      path("${datasetID}/*"), optional: true

  script:
  """
  mkdir ${datasetID}
  metabat2 \
    -i megahit/${datasetID}.contigs.fa \
    -a ${datasetID}_depth.txt \
    -o ${datasetID}/${datasetID}.individ.bin \
    -t $task.cpus \
    -m 2000 \
    -v
  """
}

process CHECKM_SINGLE {

  label 'cpus_medium'
  
  publishDir "$projectDir/results/indiv_assemblies/checkM2_output"
  
  input:
    path db
    tuple \
      val(datasetID), \
      path (metabat2_individ_outfiles, stageAs: "indiv_assembled_bins/*")

  output:
      tuple \
        val(datasetID), \
        path("$datasetID/*")
 
  script:
  """
  export HDF5_USE_FILE_LOCKING='FALSE'
  mkdir ${datasetID}
  checkm2 predict \
    --database_path $db \
    --threads 20 \
    -x fa \
    --input indiv_assembled_bins \
    --output-directory ${datasetID}
  """
}



process GET_BINS {

  label 'HQ_bins'
  publishDir "$projectDir/results/bins"
  
  input:
      path (tsv_files, stageAs: "checkM2_hq/*")
      path (individ_assembled_bins, stageAs: "bins/*")
      path (coassembled_bins, stageAs: "bins/*")
      
  output:
      path("High_quality_bins.txt")
      path("hq_bins/*")
      path("all_bins/*")
      
  script:
  """
  mkdir hq_bins
  mkdir all_bins
  
  cat checkM2_hq/*.tsv > High_quality_bins.txt
  
  cd bins
  
  for i in `cut -f 1 ../High_quality_bins.txt` ; do l=\$(readlink \$i.fa); ln -s \$l ../hq_bins ; done
  
  cd ..
  
  # a bit cray recopy the links in another folder..
  cd bins
  for i in `ls *.fa`; do l=\$(readlink \$i); ln -s \$l ../all_bins; done
  """
}





process SORT_BINS {

  label 'bins'
  publishDir "$projectDir/results/sorted_bins"
  
  input:
    tuple \
      val(datasetID), \
      path (checkm, stageAs: "checkm2_out/*")

  output:   
      path("${datasetID}_checkM2_hq.tsv")
 
  script:
  """
  awk '{if (\$2 > 90 && \$3 < 5) {print}}' checkm2_out/quality_report.tsv >  ${datasetID}_checkM2_hq.tsv
  """
}

process SORT_BINS2 {

  label 'bins'
  publishDir "$projectDir/results/sorted_bins"
  
  input:
    tuple \
      val(datasetID), \
      path (checkm, stageAs: "checkm2_out/*")
      
  output:   
      path("${datasetID}_checkM2_i_hq.tsv")
 
  script:
  """
  awk '{if (\$2 > 90 && \$3 < 5) {print}}' checkm2_out/quality_report.tsv >  ${datasetID}_checkM2_i_hq.tsv
  """
}


process DREP {

  label 'cpus_xxlarge'
  publishDir "$projectDir/results/drep"

  input:
    path hq_bins_file
    path (hq_bins, stageAs: "hq_bins/*")
    path (all_bins, stageAs: "all_bins/*")
    
  output:   
    path("dRep_output/*")

  script:
  """
  # Tweak the High quality bins file to be used by drep
  echo "genome,completeness,contamination,strain_heterogeneity" > header
  awk {'print \$1".fa,"\$2","\$3","0'} $hq_bins_file > corpus
  cat header corpus > checkM_results.csv
  
  
  dRep dereplicate -g hq_bins/*.fa \
    -comp 90 -con 5 --processors $task.cpus \
    -strW 1 -pa 0.90 -sa 0.99 \
    --S_algorithm fastANI \
    --multiround_primary_clustering \
    --greedy_secondary_clustering \
    --run_tertiary_clustering \
    --genomeInfo checkM_results.csv dRep_output

  """
}


process QUAST {

  label 'cpus_xlarge'
  
  publishDir "$projectDir/results/quast"
  
  input:
    path(coassembly, stageAs: "Megahit_coassembly/*")
    path(dereplicated_genomes, stageAs: "dRep_output/*")


  output:
    path "QUAST_replicated_MAGs/*"
    
    
    
  script:
  """
  quast.py dRep_output/dereplicated_genomes/*.fa \
    --threads $task.cpus \
    -o QUAST_replicated_MAGs
  
  metaquast.py Megahit_coassembly/Coassembly.contigs.fa \
    -t $task.cpus \
    -o QUAST_coassembly
  """


}



process GTDB_TK {

  label 'cpus_xxlarge'
  
  publishDir "$projectDir/results/GTDB"

  input:
    path(db)
    path(dereplicated_genomes, stageAs: "dRep_output/*")
 
  output:
    path("GTDBtk_output/*")
  
  
  script:
  """
  export GTDBTK_DATA_PATH=$db
  
  gtdbtk classify_wf \
         --genome_dir \
         dRep_output/dereplicated_genomes \
         -x fa \
         --out_dir GTDBtk_output \
         --cpus $task.cpus
  """
}


process PHYLOPHLAN {

 label 'cpus_xxlarge'
 
 publishDir "$projectDir/results/phylophlan"

  input:
      path(dereplicated_genomes, stageAs: "dRep_output/*")
 
  output:   
      path("Phylophlan_output/*")
  
  script:
  """
  # To avoid errors like this :
  # [e] database directory "phylophlan_databases/" is not writeable, please modify the permissions
  # We create the 'phylophlan_databases' folder before Phylophlan and change the permissions

  mkdir -p phylophlan_databases &&
  chmod a+rwx phylophlan_databases &&

  # for generating the four default configuration files
  phylophlan_write_default_configs.sh &&
  
  phylophlan -d phylophlan \
             -i dRep_output/dereplicated_genomes \
             -o Phylophlan_output \
             --db_type a \
             -f supermatrix_aa.cfg \
             --nproc $task.cpus \
             --diversity low \
             --fast \
             --verbose \
             --genome_extension fa
  """
}


process COVERM {

 label 'cpus_xlarge'
 publishDir "$projectDir/results/coverM"

  input:
    tuple \
      val(datasetID), \
      path(final_R1), \
      path(final_R2)
    path(dereplicated_genomes, stageAs: "dRep_output/*")
 
  output:
     tuple \
       val(datasetID), \
       path("${datasetID}_coverM_output.txt")
  
  script:
  """
  coverm genome -1 ${final_R1} \
                -2 ${final_R2} \
                --genome-fasta-directory dRep_output/dereplicated_genomes \
                --genome-fasta-extension fa \
                --min-covered-fraction 1 \
                --threads $task.cpus -v \
                --output-file ${datasetID}_coverM_output.txt
  """
}

process KRAKEN2 {

label 'cpus_xlarge'
label 'mem_large'
publishDir "$projectDir/results/kraken2"

input:
  path db
  tuple \
    val(datasetID), \
    path(final_R1), \
    path(final_R2)
 
output:   
  tuple \
    val(datasetID), \
    path("Kraken2_${datasetID}.report.txt")

script:
"""
kraken2 --use-names \
--threads $task.cpus \
--db $db \
--paired ${final_R1} ${final_R2} \
--report Kraken2_${datasetID}.report.txt \
--report-zero-counts > /dev/null
# It is important to redirect the large Kraken2 output to /dev/null
# Otherwise, massive info is written in .command.log
# and .command.out Nextflow files
"""
}

process KRAKEN2_MPA {

label 'cpus_xlarge'
label 'mem_large'
publishDir "$projectDir/results/kraken2_mpa"

input:
  path db
   tuple \
    val(datasetID), \
    path(final_R1), \
    path(final_R2)
 
output:   
  tuple \
    val(datasetID), \
    path("Kraken2_${datasetID}.mpa.report.txt")

script:
"""
kraken2 --use-names \
--threads $task.cpus \
--db $db \
--report Kraken2_${datasetID}.mpa.report.txt \
--use-mpa-style \
--report-zero-counts \
--paired ${final_R1} ${final_R2} > /dev/null
# It is important to redirect the large Kraken2 output to /dev/null
# Otherwise, massive info is written in .command.log
# and .command.out Nextflow files
"""
}



process COMBINE_KRAKEN2 {

label 'cpus_large'
publishDir "$projectDir/results/kraken2_summary"

input:
  path (reports, stageAs: "reports/*")

output:
  path ("Combined_Kraken2.reports.txt")
  
script:
"""
$baseDir/src/combine_mpa.py \
         -i reports/*.report.txt \
         -o Combined_Kraken2.reports.txt
"""
}


process BRACKEN {

label 'cpus_large'
publishDir "$projectDir/results/bracken"

input:
  path db
  tuple \
    val(datasetID), \
    path(kraken_files), \
    path(report_files)

output:
  tuple \
    val(datasetID), \
    path ("${datasetID}_bracken_report_species.txt")
  
script:
"""
bracken -d $db \
        -i ${report_files} \
        -r 150 -t 10 -l S \
        -o ${datasetID}_bracken_report_species.txt
"""
} 



/*
 * This process uses a custom script for producing Bracken files at various taxonomy levels
 * Author: Xavier Monger (Anthony Vicent's lab), adapted by Jean-Simon Brouard
 */
 
process BRACKEN_ALT {

label 'cpus_medium'
publishDir "$projectDir/results/bracken_smart"

input:
  path db
  path (dereplicated_genomes, stageAs: "k2_assembly_reports/*")

output:
  path ("braken/*")
  path ("bracken_abundance_files/*")

script:
"""
# produce_bracken_nf.sh call braken which is the path of the kraken2 conda env
# but it calls also 2 scripts that are located in the src directory
# We therefore specify the location of the kraken2 db and the location of these
# scripts using the built-in variable baseDir

$baseDir/src/produce_bracken_nf.sh $params.kraken2 $baseDir/src
"""
}



process DRAM_ANNOTATION {
label 'cpus_xxlarge'
label 'mem_large'

publishDir "$projectDir/results/dram/annotation"

input:
  path dram_config
  path (dereplicated_genomes, stageAs: "dRep_output/*")
  path (GTDB, stageAs: "GTDB_TK_output/*")
  

output:
  path ("DRAM_annotated_MAGs/*"), emit: DRAM_MAGs
  path 'config_infos.txt'

script:
"""
# If this file exists..
if [ -f 'gtdbtk.ar53.summary.tsv' ]; then
    tail -n +2 GTDB_TK_output/gtdbtk.ar53.summary.tsv > archae
    cat GTDB_TK_output/gtdbtk.bac120.summary.tsv archae > gtdbtk.bac120.ar53.summary.tsv
# if no archae file is present; work only with the bacterial file...
else
    mv GTDB_TK_output/gtdbtk.bac120.summary.tsv gtdbtk.bac120.ar53.summary.tsv
fi

DRAM-setup.py import_config --config_loc $dram_config

DRAM-setup.py print_config > config_infos.txt

DRAM.py annotate \
  -i 'dRep_output/dereplicated_genomes/*.fa' \
  -o DRAM_annotated_MAGs \
  --verbose \
  --config_loc $dram_config \
  --threads $task.cpus \
  --gtdb_taxonomy gtdbtk.bac120.ar53.summary.tsv
"""
}



process DRAM_DISTILLATION {

publishDir "$projectDir/results/dram/distillation"

input:
  path (annots, stageAs: "DRAM_annotated_MAGs/*")
  
output:
  path ("MAG_DRAM_distilled_summaries/*")

script:
"""
DRAM.py distill \
  -i DRAM_annotated_MAGs/annotations.tsv \
  -o MAG_DRAM_distilled_summaries \
  --trna_path DRAM_annotated_MAGs/trnas.tsv \
  --rrna_path DRAM_annotated_MAGs/rrnas.tsv \
  --genomes_per_product 575
"""
}





























