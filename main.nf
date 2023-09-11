/* 
 * This workflow is adpated from  the work of Devin Holman and Arun Kommadath (AAFC Lacombe)
 * and the work of Sara Ricci a post-doc in Renee Petri team
 *
 * 
 */

/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

log.info """\
M E T A G E N O M I C  -  Workflow - AAFC    v 0.1 
================================
genome   : $params.genome
genome basename : $params.genome_basename
reads    : $params.reads
results  : $params.results
"""

include {
  RENAME_SEQUENCES;
  QUALITY_FILTERING;
  BOWTIE2;
  OUTPUT_UNALIGNED_READS;
  KAIJU;
  KAIJU_TAX_TABLE;
  KAIJU_FULL_TAX_TABLE;
  MERGE_TAX_FILES;
  CAT_FASTQ;
  HUMANN_RUN;
  HUMANN_ABUNDANCE;
  COASSEMBLY;
  BOWTIE2_BUILD;
  BOWTIE2_MAP;
  METABAT2_BIN_COASSEMBLY;
  JGI_SUMMARIZE;
  CHECKM;
  MEGAHIT_SINGLE;
  BOWTIE2_BUILD_SINGLE;
  BOWTIE2_MAP_SINGLE;
  JGI_SUMMARIZE_SINGLE;
  METABAT2_BIN_SINGLE;
  CHECKM_SINGLE;
  SORT_BINS;
  SORT_BINS2;
  GET_BINS;
  DREP;
  GTDB_TK;
  PHYLOPHLAN;
  COVERM;
  QUAST;
  KRAKEN2;
  KRAKEN2_MPA;
  COMBINE_KRAKEN2;
  BRACKEN;
  BRACKEN_ALT;
  DRAM_ANNOTATION;
  DRAM_DISTILLATION} from './modules.nf'

  include { clean_work_files } from './utilities.nf'
/* 
 * sub workflows
 */

/*
 * sub workflow rename: This workflow rename sequences by sample id given a map_file is given.
 */
workflow rename {
    
    println "You are using the *rename* subworkflow."
    take: data
          map_file
    main:
      RENAME_SEQUENCES(data.collect(), map_file)
      
    emit:
      RENAME_SEQUENCES.out
}

/*
 * sub workflow get_reads_pairs: This workflow creates the `read_pairs_ch`
 * channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file.
 */
workflow get_reads_pairs {
    
    println "You are using the *get_reads_pairs* subworkflow."
    take: data
    
    main:
       data.flatten()
       .filter( ~/^.*fastq.gz/ )
       .map { file ->
          def key_part1 = file.name.toString().tokenize('_').get(0)
          key = key_part1
          return tuple(key, file)
       }
       .groupTuple(size:2)
       .flatten()
       .collate ( 3 )
       .set { read_pairs_ch }
      
    emit:
      read_pairs_ch
}


/*
 * sub workflow get_folder: This workflow creates a `folder_ch`
 * channel that emits tuples containing 2 elements:
 * the sample ID, and the path containing a folder with stuff.
 */
workflow get_folder {
    
    println "You are using the *get_folder* subworkflow."
    take: data
    
    main:
       data.flatten()
       .map { dir ->
          def key_part1 = dir.name.toString()
          key = key_part1
          return tuple(key, dir)
       }
       .view()       
       .set { folder_ch }
    emit:
      folder_ch
}


/* 
 * main pipeline logic
 */
workflow {

  if( params.use_prepared_reads ) {

    // Using the get_reads_pairs workflow
    prepared_reads_ch = get_reads_pairs( Channel.fromPath( params.prepared_reads) )


  } else {

    // using the rename workflow with 2 inputs
    rename( Channel.fromPath( params.reads), params.map_file )

    // using the get_reads_pairs workflow
    get_reads_pairs(rename.out)
    
    // PART 1: Data preparation
    QUALITY_FILTERING(get_reads_pairs.out)
    BOWTIE2(params.genome, params.genome_basename,
            QUALITY_FILTERING.out)
            
    
    prepared_reads_ch = OUTPUT_UNALIGNED_READS(BOWTIE2.out)
    
    // This signal will triggered the clean files process
    QUALITY_FILTERING.out.join(OUTPUT_UNALIGNED_READS.out).join(BOWTIE2.out)
      .collect()
      .flatten()
      //.filter{ it =~ /.*bam$/ }
      .filter{ it =~ /.*_trimmed.*/ || it =~ /.*bam$/ || it =~ /.*_unpaired.*/ }
      .view()
      .set{ cleanable_bams_ch }
    clean_sorted_bams(cleanable_bams_ch)
  }


    KAIJU(params.kaiju_db, prepared_reads_ch)
    ch_kaiju = KAIJU_TAX_TABLE(params.kaiju_db,KAIJU.out)
    KAIJU_FULL_TAX_TABLE(params.kaiju_db,KAIJU.out)   

    ch_kaiju
      .flatten()
      .filter ( Path ) // To get rid of datasetID values    
      .collect()     
      //.view()
      .set { ch_kaiju_tsv }    

    MERGE_TAX_FILES(ch_kaiju_tsv)
    CAT_FASTQ(prepared_reads_ch)
    HUMANN_RUN(params.chocophlan_db, params.metaphlan_db, params.uniref_db, CAT_FASTQ.out)
    HUMANN_ABUNDANCE(HUMANN_RUN.out.flatten().filter ( Path ).collect())   
    COASSEMBLY(prepared_reads_ch.flatten().filter ( ~/^.*R1.fastq.gz/ ).collect(),
               prepared_reads_ch.flatten().filter ( ~/^.*R2.fastq.gz/ ).collect())
    BOWTIE2_BUILD(COASSEMBLY.out)
    BOWTIE2_MAP(BOWTIE2_BUILD.out,prepared_reads_ch)    
    JGI_SUMMARIZE(BOWTIE2_MAP.out)
    METABAT2_BIN_COASSEMBLY(COASSEMBLY.out,JGI_SUMMARIZE.out)
    CHECKM(params.checkm2_db, METABAT2_BIN_COASSEMBLY.out)
    

    if( params.use_megahit_individual_assemblies ) {

    // Using the get_folders workflow
    indiv_assemblies_ch = get_folder( Channel.fromPath( params.indiv_assemblies, type: 'dir') )

    } else {
    // PART 2 : Individual assemblies
    indiv_assemblies_ch = MEGAHIT_SINGLE(prepared_reads_ch)
    
    }

    
    BOWTIE2_BUILD_SINGLE(indiv_assemblies_ch)
    BOWTIE2_MAP_SINGLE(BOWTIE2_BUILD_SINGLE.out.join(prepared_reads_ch))
    JGI_SUMMARIZE_SINGLE(BOWTIE2_MAP_SINGLE.out)
    ch_meta = METABAT2_BIN_SINGLE(JGI_SUMMARIZE_SINGLE.out.join(indiv_assemblies_ch))
    
    
    CHECKM_SINGLE(params.checkm2_db, METABAT2_BIN_SINGLE.out)
    SORT_BINS(CHECKM.out)
    SORT_BINS2(CHECKM_SINGLE.out)
    GET_BINS(SORT_BINS.out.concat(SORT_BINS2.out).collect(),
             METABAT2_BIN_SINGLE.out.flatten().filter ( Path ).collect(),
             METABAT2_BIN_COASSEMBLY.out.flatten().filter ( Path ).collect())

    DREP(GET_BINS.out)
    QUAST(COASSEMBLY.out, DREP.out)
    GTDB_TK(params.gtdb_db, DREP.out)
    PHYLOPHLAN(DREP.out)
    COVERM(prepared_reads_ch,DREP.out)
    
    KRAKEN2(params.kraken2, prepared_reads_ch)
    KRAKEN2_MPA(params.kraken2, prepared_reads_ch)
    COMBINE_KRAKEN2(KRAKEN2_MPA.out.flatten().filter ( Path ).collect())
    BRACKEN_ALT(params.kraken2, KRAKEN2.out.flatten().filter ( Path ).collect())

    DRAM_ANNOTATION(params.dram_config, DREP.out, GTDB_TK.out)
    DRAM_DISTILLATION(DRAM_ANNOTATION.out.DRAM_MAGs)
}



