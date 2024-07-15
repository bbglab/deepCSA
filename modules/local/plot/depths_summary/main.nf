process PLOT_DEPTHS {
    tag "$meta.id"
    label 'process_single'
    label 'time_low'
    label 'process_high_memory'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seaborn:0.12.2_cv1' :
        'biocontainers/seaborn:0.12.2_cv1' }"

    input:
    tuple val(meta) , path(depth)
    tuple val(meta2), path(panel)

    output:
    tuple val(meta), path("*.pdf")          , emit: plots
    tuple val(meta), path("*depth*.tsv")    , emit: depths
    path  "versions.yml"                    , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    plot_depths.py \\
                ${prefix} \\
                ${depth} \\
                ${panel} \\
                ${panel_version};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    touch ${prefix}.${panel_version}.depths_info.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
