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
  QUALITY_FILTERING;
  BOWTIE2;
  SAM2BAM;
  SORTBAM;
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
  CHECKM_SINGLE} from './modules.nf'

    

/* 
 * main pipeline logic
 */
workflow {

    // PART 1: Preparation of the raw reads channel
    // Create the `read_pairs_ch` channel that emits tuples containing three elements:
    // the pair ID, the first read-pair file and the second read-pair file
    
    Channel
        .fromFilePairs( params.reads, flat:true)
        .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
        //.view()
        .set { read_pairs_ch } 

    // PART 1: Data preparation
    QUALITY_FILTERING(read_pairs_ch)
    BOWTIE2(params.genome, params.genome_basename,
            QUALITY_FILTERING.out)
    SAM2BAM(BOWTIE2.out)
    SORTBAM(SAM2BAM.out)
    OUTPUT_UNALIGNED_READS(SORTBAM.out)
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
    
    
    CHECKM(METABAT2_BIN_COASSEMBLY.out)
    
    // PART 2 : Individual assemblies
    MEGAHIT_SINGLE(OUTPUT_UNALIGNED_READS.out)
    BOWTIE2_BUILD_SINGLE(MEGAHIT_SINGLE.out)
    BOWTIE2_MAP_SINGLE(BOWTIE2_BUILD_SINGLE.out.join(OUTPUT_UNALIGNED_READS.out))
    SORTSAM_SINGLE(BOWTIE2_MAP_SINGLE.out)
    JGI_SUMMARIZE_SINGLE(SORTSAM_SINGLE.out)
    METABAT2_BIN_SINGLE(JGI_SUMMARIZE_SINGLE.out.join(MEGAHIT_SINGLE.out))
    CHECKM_SINGLE(METABAT2_BIN_SINGLE.out)


}


/*
    // hack to test the CHECKM_SINGLE module
     Channel.fromPath(params.pseudobins)
    .map { file ->
        def key_part1 = file.name.toString().split('.individ.bin.1.fa')
        key = key_part1
        return tuple(key, file)
     }
    .flatten()
    .collate (2)
    .view()
    .set{ pseudobins_ch }

*/




