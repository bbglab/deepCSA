process HOTSPOTS_SELECTION {
    tag "$meta.id"

    label 'deepcsa_core'

    input:
    tuple val(meta), path(comparisons)
    tuple val(meta2), path(annotated_panel_richer)
    path(hotspots_file)

    output:
    tuple val(meta), path("*.hotspots_selection.tsv.gz") , emit: selections
    path "versions.yml"                                  , topic: versions

    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def comparison_args = comparisons.collect { "--comparisons ${it}" }.join(' ')
    """
    compute_hotspots_selection.py \\
        ${comparison_args} \\
        --panel-file ${annotated_panel_richer} \\
        --hotspots-file ${hotspots_file} \\
        --output-prefix ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.site.hotspots_selection.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
