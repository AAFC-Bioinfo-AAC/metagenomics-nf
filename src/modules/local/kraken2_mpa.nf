process KRAKEN2_MPA {

    publishDir "$params.results/kraken2_mpa"

    input:
    tuple \
        val(datasetID), \
        path(report)
    
    output:   
    tuple \
        val(datasetID), \
        path("Kraken2_${datasetID}.mpa.report.txt")

    script:
    """
    $baseDir/src/kreport2mpa.py \
    -r $report \
    -o Kraken2_${datasetID}.mpa.report.txt
    """

}