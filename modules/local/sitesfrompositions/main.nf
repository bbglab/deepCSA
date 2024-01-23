process SITESFROMPOSITIONS {
    // TODO
    // add proper tags

    tag "${meta.id}"
    // label 'cpu_single'
    // label 'time_low'
    // label 'process_low_memory'

    container "docker.io/ferriolcalvet/bgreference"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //         'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
    //         'biocontainers/pandas:1.5.2' }"


    input:
    tuple val(meta), path(depths)

    output:
    path("*.tsv")        , emit: annotated_panel_reg
    path  "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "TargetRegions"
    """
    zcat ${depths} | cut -f1,2  > captured_positions.tsv;
    sites_table_from_positions.py \\
                    captured_positions.tsv \\
                    captured_positions.sites4VEP.tsv;
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "TargetRegions"
    """
    touch captured_positions.sites4VEP.tsv;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
