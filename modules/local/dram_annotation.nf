process DRAM_ANNOTATION {

    label 'cpus_xxlarge'
    label 'mem_large'

    publishDir "$params.results/dram/annotation"

    input:
    path config
    path (dereplicated_genomes, stageAs: "dRep_output/*")
    path (GTDB, stageAs: "GTDB_TK_output/*")
    

    output:
    path ("DRAM_annotated_MAGs/*"), emit: DRAM_MAGs

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

    export DRAM_CONFIG_LOCATION=$config

    DRAM-setup.py print_config > config_infos.txt
    
    DRAM.py annotate \
    -i 'dRep_output/*.fa' \
    -o DRAM_annotated_MAGs \
    --verbose \
    --config_loc $config \
    --threads $task.cpus \
    --gtdb_taxonomy gtdbtk.bac120.ar53.summary.tsv
    """
}
