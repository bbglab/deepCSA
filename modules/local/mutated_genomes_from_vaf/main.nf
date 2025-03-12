process MUTATED_GENOMES_FROM_VAF {
    tag "$meta.id"
    label 'process_single'
    label 'memory_medium'

    container 'docker.io/ferranmuinos/test_mutated_genomes'

    input:
    tuple val(meta) , path(mutations), path(omegas)

    output:
    tuple val(meta), path("*.covered_genomes_summary.tsv") , emit: mutated_gen_sample
    path  "versions.yml"                                   , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mutgenomes_driver_priority.py \\
                --sample ${prefix} \\
                --somatic-mutations-file ${mutations} \\
                --omega-file ${omegas} ;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    def panel_version = task.ext.panel_version ?: "${meta.id}"
    """
    touch ${prefix}.${panel_version}.mutrates.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
