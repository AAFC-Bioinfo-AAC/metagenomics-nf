process DRAM_ANNOTATION {

    label 'cpus_xxlarge'
    label 'mem_large'

    publishDir "$params.results/dram/annotation"

    input:
    path dram_config
    path (dereplicated_genomes, stageAs: "dRep_output/*")
    path (GTDB, stageAs: "GTDB_TK_output/*")
    

    output:
    path ("DRAM_annotated_MAGs/*"), emit: DRAM_MAGs
    path 'config_infos.txt'

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

    DRAM-setup.py import_config --config_loc $dram_config

    DRAM-setup.py print_config > config_infos.txt

    DRAM.py annotate \
    -i 'dRep_output/dereplicated_genomes/*.fa' \
    -o DRAM_annotated_MAGs \
    --verbose \
    --config_loc $dram_config \
    --threads $task.cpus \
    --gtdb_taxonomy gtdbtk.bac120.ar53.summary.tsv
    """
    
}