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

/*
 * Define the default parameters
 */ 
params.genome     = "$baseDir/data/genome.fa"
params.reads      = "$baseDir/data/raw_reads/*.fastq"
params.results    = "results"


log.info """\
M E T A G E N O M I C  -  Workflow - AAFC    v 0.1 
================================
genome   : $params.genome
reads    : $params.reads
results  : $params.results
"""

include { 
  QUALITY_FILTERING;
  CALLING_VARIANTS_PER_SAMPLE} from './modules.nf'


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

}
