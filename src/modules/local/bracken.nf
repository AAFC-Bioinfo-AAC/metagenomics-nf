/*
 * This process uses a custom script for producing Bracken files at various taxonomy levels
 * Author: Xavier Monger (Anthony Vicent's lab), adapted by Jean-Simon Brouard
 */

process BRACKEN {

    label 'cpus_medium'
    publishDir "$params.results/bracken_smart"

    input:
    path db
    path (kraken2_reports, stageAs: "k2_assembly_reports/*")

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