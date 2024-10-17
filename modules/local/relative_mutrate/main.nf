process RELATIVE_MUTRATE {
    tag "$meta.id"
    label 'process_single'

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
    //     'biocontainers/pandas:1.5.2' }"
    container 'docker.io/ferriolcalvet/bgreference'


    input:
    tuple val(meta), path(mutation_rates)

    output:
    tuple val(meta), path("*.rel_mutrates.tsv") , emit: mutrate
    path  "versions.yml"                        , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    omega_compute_relative_mutrate.py \\
                --mutrates ${mutations} \\
                --output ${prefix}.rel_mutrates.tsv;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    touch ${prefix}.rel_mutrates.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
