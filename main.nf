/* 
 * This workflow was adpated by Jean-Simon Brouard
 * from the work done by Devin Holman, Arun Kommadath (AAFC Lacombe) and Sara Ricci.
 * Last update : 2023/11/08
 * 
*/

/* 
 * Enable DSL 2 syntax
*/
nextflow.enable.dsl = 2

log.info """\
M E T A G E N O M I C  -  Workflow - AAFC    v 1.1 
================================
genome   : $params.genome
genome basename : $params.genome_basename
reads    : $params.reads
results  : $params.results
"""

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO: Use nf-core modules where possible for easier updating
// TODO: Consolidate modules where possible (e.g. Bowtie2)
include { BOWTIE2_BUILD_SINGLE            } from './modules/local/bowtie2_build_single'
include { BOWTIE2_BUILD                   } from './modules/local/bowtie2_build'
include { BOWTIE2_MAP_SINGLE              } from './modules/local/bowtie2_map_single'
include { BOWTIE2_MAP                     } from './modules/local/bowtie2_map'
include { BOWTIE2                         } from './modules/local/bowtie2'
include { BRACKEN                         } from './modules/local/bracken'
include { CAT_FASTQ                       } from './modules/local/cat_fastq'
include { CHECKM_SINGLE                   } from './modules/local/checkm_single'
include { CHECKM                          } from './modules/local/checkm'
include { COASSEMBLY                      } from './modules/local/coassembly'
include { COMBINE_KRAKEN2                 } from './modules/local/combine_kraken2'
include { COVERM                          } from './modules/local/coverm'
include { DRAM_ANNOTATION                 } from './modules/local/dram_annotation'
include { DRAM_DISTILLATION               } from './modules/local/dram_distillation'
include { DREP                            } from './modules/local/drep'
include { GET_BINS                        } from './modules/local/get_bins'
include { GET_BINS2                       } from './modules/local/get_bins2'
include { GTDB_TK                         } from './modules/local/gtdb_tk'
include { HUMANN_ABUNDANCE                } from './modules/local/humann_abundance'
include { HUMANN_RUN                      } from './modules/local/humann_run'
include { JGI_SUMMARIZE_SINGLE            } from './modules/local/jgi_summarize_single'
include { JGI_SUMMARIZE                   } from './modules/local/jgi_summarize'
include { KAIJU_FULL_TAX_TABLE            } from './modules/local/kaiju_full_tax_table'
include { KAIJU_TAX_TABLE                 } from './modules/local/kaiju_tax_table'
include { KAIJU                           } from './modules/local/kaiju'
include { KRAKEN2_MPA                     } from './modules/local/kraken2_mpa'
include { KRAKEN2                         } from './modules/local/kraken2'
include { MEGAHIT_SINGLE                  } from './modules/local/megahit_single'
include { MERGE_FULL_TAX_FILES            } from './modules/local/merge_full_tax_files'
include { MERGE_TAX_FILES                 } from './modules/local/merge_tax_files'
include { METABAT2_BIN_COASSEMBLY         } from './modules/local/metabat2_bin_coassembly'
include { METABAT2_BIN_SINGLE             } from './modules/local/metabat2_bin_single'
include { METAQUAST                       } from './modules/local/metaquast'
include { OUTPUT_UNALIGNED_READS          } from './modules/local/output_unaligned_reads'
include { PHYLOPHLAN                      } from './modules/local/phylophlan'
include { QUALITY_FILTERING               } from './modules/local/quality_filtering'
include { QUAST                           } from './modules/local/quast'
include { SORT_BINS                       } from './modules/local/sort_bins'
include { SORT_BINS2                      } from './modules/local/sort_bins2'
// TODO: Make the GET_FOLDER, GET_READS_PAIRS, and RENAME subworkflows normal modules?
include { GET_FOLDER                      } from './subworkflows/local/get_folder'
include { GET_READS_PAIRS                 } from './subworkflows/local/get_reads_pairs'
include { RENAME                          } from './subworkflows/local/rename'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT UTILITIES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO: Make this a normal module?
include { clean_work_files as clean_sorted_bams } from './utilities'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN PIPELINE LOGIC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

  if( params.use_prepared_reads ) {

    // Using the GET_READS_PAIRS workflow
    prepared_reads_ch = GET_READS_PAIRS( Channel.fromPath( params.prepared_reads) )


  } else {

    // Using the RENAME subworkflow with 2 inputs
    RENAME( Channel.fromPath( params.reads), params.map_file )

    // Using the GET_READS_PAIRS workflow
    GET_READS_PAIRS(RENAME.out)

    // PART 1: Data preparation
    QUALITY_FILTERING(GET_READS_PAIRS.out)
    BOWTIE2(params.genome, params.genome_basename,
            QUALITY_FILTERING.out)

    prepared_reads_ch = OUTPUT_UNALIGNED_READS(BOWTIE2.out)
    
    // This signal will trigger the clean files process
    QUALITY_FILTERING.out.join(OUTPUT_UNALIGNED_READS.out).join(BOWTIE2.out)
      .collect()
      .flatten()
      //.filter{ it =~ /.*bam$/ }
      .filter{ it =~ /.*_trimmed.*/ || it =~ /.*bam$/ || it =~ /.*_unpaired.*/ }
      //.view()
      .set{ cleanable_bams_ch }
    clean_sorted_bams(cleanable_bams_ch)
  }

  if (!params.skip_kaiju ) {
    println "*You do Kaiju*"
    KAIJU(params.kaiju_db, prepared_reads_ch)
    //ch_kaiju = KAIJU_TAX_TABLE(params.kaiju_db, KAIJU.out)
    ch_kaiju_full = KAIJU_FULL_TAX_TABLE(params.kaiju_db, KAIJU.out)   
    //ch_kaiju
    //  .flatten()
    //  .filter ( Path ) // To get rid of datasetID values    
    //  .collect()     
    //  .set { ch_kaiju_tsv }
    ch_kaiju_full
      .flatten()
      .filter ( Path ) // To get rid of datasetID values    
      .collect()     
      .set { ch_kaiju_full_tsv }
      //MERGE_TAX_FILES(ch_kaiju_tsv)
      MERGE_FULL_TAX_FILES(ch_kaiju_full_tsv)

    } else {
      println "*You skip Kaiju...*"
  }

  if (!params.skip_kraken ) {
    println "*You do Kraken2*"
    KRAKEN2(params.kraken2_confidence_threshold, params.kraken2, prepared_reads_ch)
    KRAKEN2_MPA(KRAKEN2.out)
    COMBINE_KRAKEN2(KRAKEN2_MPA.out.flatten().filter ( Path ).collect())
    BRACKEN(params.kraken2, KRAKEN2.out.flatten().filter ( Path ).collect())
  } else {
    println "*You skip Kraken...*"
  }  

  if (!params.skip_humann ) {
    println "*You do Humann*"
    CAT_FASTQ(prepared_reads_ch)
    HUMANN_RUN(params.chocophlan_db, params.metaphlan_db, params.uniref_db, CAT_FASTQ.out)
    HUMANN_ABUNDANCE(HUMANN_RUN.out.flatten().filter ( Path ).collect())   
  } else {
    println "*You skip Humann...*"
  }   

  // Individual assemblies
  if( params.use_megahit_individual_assemblies ) {
    // Using the GET_FOLDER subworkflow
    indiv_assemblies_ch = GET_FOLDER( Channel.fromPath( params.indiv_assemblies, type: 'dir') )
  } else {
    indiv_assemblies_ch = MEGAHIT_SINGLE(prepared_reads_ch)
  }
  
  if (!params.skip_coassembly ) {
    println "*You do both individual and coassembly steps*"

    /*
    * Create the `type_ch` channel using the 'type' column from the metadata
    * Example :
    *[66b, rumen]
    *[67b, milk]
    *[68b, milk]
    */
    Channel
      .fromPath(params.coassembly_file)
      .splitCsv(header: true, sep: "\t")
      .map{ row-> 
            def key_part1 = row.sample_read.toString().tokenize('_').get(0)
            return tuple(key_part1, "${row.coassembly_group}")
          }
      .distinct()
      .set { type_ch }

    /*
    * Create the `reads_plus_ch` channel
    * First we join the information of the sample type to the reads_ch
    * Both channels have the sampleID as primary key
    * Then we move this information at the beginning of the tuple so it can be
    * used as a new key
    * After a groupTuple allows to regroup all reads belonging to the same type
    * Example :
    * [rumen, [.../renamed/66b_R1.fastq.gz], [.../renamed/66b_R2.fastq.gz], [66b]]
    * [milk, [.../renamed/67b_R1.fastq.gz, .../68b_R1.fastq.gz], [.../renamed/67b_R2.fastq.gz, .../renamed/68b_R2.fastq.gz], [67b, 68b]]
    */
    prepared_reads_ch.join(type_ch)
                     .map { row-> 
                              def a = row[0]
                              def b = row[1]
                              def c = row[2]
                              def d = row[3]
                            return tuple(d,b,c,a)
                          }
                     .groupTuple()
                     .set {reads_plus_ch}
    
    /*
    * Create the `prepared_reads_like_ch` channel
    * First we join the information of the sample type to the reads_ch
    * Both channels have the sampleID as primary key
    * Then we move this information at the beginning of the tuple so it can be
    * used as a new key for joining with BOWTIE2_BUILD
    */
    prepared_reads_ch.join(type_ch)
                     .map { row-> 
                              def a = row[0]
                              def b = row[1]
                              def c = row[2]
                              def d = row[3]
                            return tuple(d,b,c,a)
                          }
                     .set {prepared_reads_like_ch}

    //prepared_reads_like_ch.view()

    COASSEMBLY(reads_plus_ch)
    BOWTIE2_BUILD(COASSEMBLY.out)

    /* Prepare the prepared_reads_and_index_ch channel containg the prepared reads
    *  and their associated co-assembly genome index. We need to use the combine operator (not join!) e.g. :
    * [rumen_7, .../C618-d7-WRC_R1.fastq.gz, .../C618-d7-WRC_R2.fastq.gz, C618-d7-WRC, [.../rumen_7/coassembly.1.bt2, ...  .../rumen_7/coassembly.rev.2.bt2]]
    * [rumen_2, .../C830-d2-WRC_R1.fastq.gz, .../C830-d2-WRC_R2.fastq.gz, C830-d2-WRC, [.../rumen_2/coassembly.1.bt2, ...  .../rumen_2/coassembly.rev.2.bt2]]
    * [milk_2,  .../C844-d2-M_R1.fastq.gz,   .../C844-d2-M_R2.fastq.gz,   C844-d2-M,   [.../milk_2/coassembly.1.bt2,  ...  .../milk_2/coassembly.rev.2.bt2]]
    * [milk_2,  .../C255-d2-M_R1.fastq.gz,   .../C255-d2-M_R2.fastq.gz,   C255-d2-M,   [.../milk_2/coassembly.1.bt2,  ...  .../milk_2/coassembly.rev.2.bt2]]
    * etc
    */
    prepared_reads_and_index_ch = prepared_reads_like_ch.combine(BOWTIE2_BUILD.out, by: 0)

    BOWTIE2_MAP(prepared_reads_and_index_ch)        
    JGI_SUMMARIZE(BOWTIE2_MAP.out)
    METABAT2_BIN_COASSEMBLY(COASSEMBLY.out.combine(JGI_SUMMARIZE.out, by: 0))
    CHECKM(params.checkm2_db, METABAT2_BIN_COASSEMBLY.out)
  
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
    QUAST(DREP.out)
    METAQUAST(COASSEMBLY.out)
    GTDB_TK(params.gtdb_db, DREP.out)
    PHYLOPHLAN(params.phylophlan_db, DREP.out)
    COVERM(prepared_reads_ch,DREP.out)
    DRAM_ANNOTATION(params.dram_config, DREP.out, GTDB_TK.out)
    DRAM_DISTILLATION(DRAM_ANNOTATION.out.DRAM_MAGs)

    } else {

    println "*You skip all coassembly steps...*"
    BOWTIE2_BUILD_SINGLE(indiv_assemblies_ch)
    BOWTIE2_MAP_SINGLE(BOWTIE2_BUILD_SINGLE.out.join(prepared_reads_ch))
    JGI_SUMMARIZE_SINGLE(BOWTIE2_MAP_SINGLE.out)
    ch_meta = METABAT2_BIN_SINGLE(JGI_SUMMARIZE_SINGLE.out.join(indiv_assemblies_ch))
    CHECKM_SINGLE(params.checkm2_db, METABAT2_BIN_SINGLE.out)
    SORT_BINS2(CHECKM_SINGLE.out)
    GET_BINS2(SORT_BINS2.out.collect(),
             METABAT2_BIN_SINGLE.out.flatten().filter ( Path ).collect())
    DREP(GET_BINS2.out)
    QUAST(DREP.out)
    GTDB_TK(params.gtdb_db, DREP.out)
    PHYLOPHLAN(params.phylophlan_db, DREP.out)
    COVERM(prepared_reads_ch,DREP.out)
    DRAM_ANNOTATION(params.dram_config, DREP.out, GTDB_TK.out)
    DRAM_DISTILLATION(params.dram_config,DRAM_ANNOTATION.out.DRAM_MAGs)

  }   
}




