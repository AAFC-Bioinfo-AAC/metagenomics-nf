process DRAM_DISTILLATION {

    publishDir "$params.results/dram/distillation"

    input:
    path (annots, stageAs: "DRAM_annotated_MAGs/*")
    
    output:
    path ("MAG_DRAM_distilled_summaries/*")

    script:
    """
    DRAM.py distill \
    -i DRAM_annotated_MAGs/annotations.tsv \
    -o MAG_DRAM_distilled_summaries \
    --trna_path DRAM_annotated_MAGs/trnas.tsv \
    --rrna_path DRAM_annotated_MAGs/rrnas.tsv \
    --genomes_per_product 575
    """
    
}