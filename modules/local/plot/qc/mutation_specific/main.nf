process PLOT_MUTATION_SPECIFIC {

    tag "${meta.id}"
    label 'process_low'

    label 'deepcsa_core'

    input:
    tuple val(meta), path (all_mutations)

    output:
    path("**.pdf")      , optional : true, emit: plots
    path("**.tsv")      , optional : true,  emit: tables
    path "versions.yml" , topic: versions

    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def max_n = task.ext.max_n ?: 10
    """
    plot_qc_mutations_vaf.py \\
                --sample_name ${prefix} \\
                --maf_file ${all_mutations} \\
                --output_prefix ${prefix} \\
                --max_n ${max_n}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
