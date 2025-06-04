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


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def recode = task.ext.recode_list ? "--recoded-genes ${task.ext.recode_list}" : ""
    """
    mutgenomes_driver_priority.py \\
                --sample ${prefix} \\
                --somatic-mutations-file ${mutations} \\
                --omega-file ${omegas} \\
                ${recode} ;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    """
    touch ${prefix}.covered_genomes_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
