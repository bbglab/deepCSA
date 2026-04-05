process ANNOTATE_DEPTHS {
    tag "${meta.id}"
    label 'process_low'
    label 'time_low'

    label 'deepcsa_core'
    
    input:
    tuple val(meta) , path(depths)
    tuple val(meta2), path(panel_all)
    path (json_groups)
    path (mask_matrix)
    path (input_csv)

    output:
    path("*.depths.annotated.tsv.gz")                           , emit: annotated_depths
    tuple val(meta), path("all_samples_indv.depths.tsv.gz")     , emit: all_samples_depths
    path "versions.yml"                                         , topic: versions


    script:
    """
    cut -f 1,2,9 ${panel_all} | uniq > ${panel_all}.contexts
    merge_annotation_depths.py \\
        --annotation ${panel_all}.contexts \\
        --depths ${depths} \\
        --json_file ${json_groups} \\
        --mask-matrix ${mask_matrix} \\
        --input-csv ${input_csv}
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch all_samples.cohort.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
