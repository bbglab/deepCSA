process FILTER_CAPTURED_PANEL {
    tag "$meta.id"
    label 'process_single'
    label 'process_medium_high_memory'

    container "community.wave.seqera.io/library/bedtools_pybedtools_pandas_pip_pruned:78080da05d53636d"

    input:
    tuple val(meta), path(complete_annotated_panel)
    tuple val(meta), path(flagged_bed)

    output:
    tuple val(meta), path("captured_panel.tab.filtered.tsv")                      , emit: filtered_annotated_panel
    path("captured_panel.tab.removed_variants.tsv")              , emit: removed_variants

    script:
    """
    # Extract header from panel and save it
    head -n 1 ${complete_annotated_panel} > captured_panel.tab.filtered.tsv
    
    # Create header for removed variants file with additional columns
    head -n 1 ${complete_annotated_panel} | awk '{print \$0 "\\tBED_CHROM\\tBED_START\\tBED_END\\tFILTER"}' > captured_panel.tab.removed_variants.tsv

    # Create temporary bed file from complete annotated panel and filter out flagged regions
    tail -n +2 ${complete_annotated_panel} | \
    awk 'BEGIN {FS=OFS="\\t"} {print \$1, \$2-1, \$2, \$0}' > temp_captured_panel.bed

    bedtools intersect -a temp_captured_panel.bed -b ${flagged_bed} -v | cut -f4- >> captured_panel.tab.filtered.tsv

    # Create a separate file with removed variants
    bedtools intersect -a temp_captured_panel.bed -b ${flagged_bed} -wa -wb | cut -f4- >> captured_panel.tab.removed_variants.tsv

    rm temp_captured_panel.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch captured_panel.tab.filtered.tsv
    touch captured_panel.tab.removed_variants.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}