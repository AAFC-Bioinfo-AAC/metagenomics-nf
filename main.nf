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
  MERGE_TAX_FILES} from './modules.nf'


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
}
