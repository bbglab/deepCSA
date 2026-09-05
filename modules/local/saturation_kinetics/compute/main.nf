process COMPUTE_SATURATION_KINETICS {

    tag "$meta.id"
    label 'cpu_low'

    label 'deepcsa_core'

    input:
    tuple val(meta) , path(mutations), path(depths)
    tuple val(meta1), path(captured_panel_rich)
    tuple val(meta2), path(expanded_panel, stageAs: "expanded_panel.tsv")
    

    output:
    tuple val(meta), path("**.tsv"), emit: table
    path "versions.yml"            , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    // --residue
    // --impact
    """
    discovery.py \\
                    --somatic-mutations-file ${mutations} \\
                    --vep-file ${captured_panel_rich} \\
                    --consensus-panel-file ${expanded_panel} \\
                    --depths-file ${depths} \\
                    --group-name ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
