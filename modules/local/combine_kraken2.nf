process COMBINE_KRAKEN2 {

    label 'cpus_large'
    publishDir "$params.results/kraken2_summary"

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