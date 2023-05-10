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
  SORTSAM;
  CHECKM;
  MEGAHIT_SINGLE;
  BOWTIE2_BUILD_SINGLE;
  BOWTIE2_MAP_SINGLE;
  SORTSAM_SINGLE;
  JGI_SUMMARIZE_SINGLE;
  METABAT2_BIN_SINGLE;
  CHECKM_SINGLE;
  SORT_BINS;
  SORT_BINS2;
  GET_BINS;
  DREP;GTDB_TK;PHYLOPHLAN;COVERM;
  QUAST;KRAKEN2;COMBINE_KRAKEN2;BRACKEN;DRAM_PREPARE_DB;DRAM_ANNOTATION} from './modules.nf'

/* 
 * sub workflows
 */



/*
 * rename: This workflow rename sequences by sample id given a map_file is given.
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
 * get_reads_pairs: This workflow creates the `read_pairs_ch` channel that emits tuples containing three elements:
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
 * main pipeline logic
 */
workflow {

    // using the rename workflow with 2 inputs
    rename( Channel.fromPath( params.reads), params.map_file )

    // using the get_reads_pairs workflow
    get_reads_pairs(rename.out)
    
  
    // PART 1: Data preparation
    QUALITY_FILTERING(get_reads_pairs.out)
    BOWTIE2(params.genome, params.genome_basename,
            QUALITY_FILTERING.out)
    OUTPUT_UNALIGNED_READS(BOWTIE2.out)
    KAIJU(params.kaiju_db, OUTPUT_UNALIGNED_READS.out)
    ch_kaiju = KAIJU_TAX_TABLE(params.kaiju_db,KAIJU.out)
    KAIJU_FULL_TAX_TABLE(params.kaiju_db,KAIJU.out)   

    ch_kaiju
      .flatten()
      .filter ( Path ) // To get rid of datasetID values    
      .collect()     
      //.view()
      .set { ch_kaiju_tsv }    

    MERGE_TAX_FILES(ch_kaiju_tsv)
    CAT_FASTQ(OUTPUT_UNALIGNED_READS.out)
    HUMANN_RUN(params.chocophlan_db, params.metaphlan_db, params.uniref_db, CAT_FASTQ.out)
    HUMANN_ABUNDANCE(HUMANN_RUN.out.flatten().filter ( Path ).collect())   
    COASSEMBLY(OUTPUT_UNALIGNED_READS.out.flatten().filter ( ~/^.*R1.fastq.gz/ ).collect(),
               OUTPUT_UNALIGNED_READS.out.flatten().filter ( ~/^.*R2.fastq.gz/ ).collect())
    BOWTIE2_BUILD(COASSEMBLY.out)
    BOWTIE2_MAP(BOWTIE2_BUILD.out,OUTPUT_UNALIGNED_READS.out)
    SORTSAM(BOWTIE2_MAP.out)
    
    JGI_SUMMARIZE(SORTSAM.out)
    METABAT2_BIN_COASSEMBLY(COASSEMBLY.out,JGI_SUMMARIZE.out)
    CHECKM(params.checkm2_db, METABAT2_BIN_COASSEMBLY.out)
    
    // PART 2 : Individual assemblies
    MEGAHIT_SINGLE(OUTPUT_UNALIGNED_READS.out)
    BOWTIE2_BUILD_SINGLE(MEGAHIT_SINGLE.out)
    BOWTIE2_MAP_SINGLE(BOWTIE2_BUILD_SINGLE.out.join(OUTPUT_UNALIGNED_READS.out))
    SORTSAM_SINGLE(BOWTIE2_MAP_SINGLE.out)
    JGI_SUMMARIZE_SINGLE(SORTSAM_SINGLE.out)
    ch_meta = METABAT2_BIN_SINGLE(JGI_SUMMARIZE_SINGLE.out.join(MEGAHIT_SINGLE.out))
    
    
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
    COVERM(OUTPUT_UNALIGNED_READS.out,DREP.out)
    
    KRAKEN2(params.kraken2, OUTPUT_UNALIGNED_READS.out)
    
    COMBINE_KRAKEN2(KRAKEN2.out.flatten().filter ( Path ).collect())
    BRACKEN(params.kraken2, KRAKEN2.out)
    
    // There is an alrady set-up database on the biocluster
    //DRAM_PREPARE_DB(params.gene_ko_link_loc, params.kegg_loc, params.viral_loc)
    DRAM_ANNOTATION(DREP.out, GTDB_TK.out)
    
    
    
}



