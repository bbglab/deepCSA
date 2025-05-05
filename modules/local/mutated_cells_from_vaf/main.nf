process MUTATED_CELLS_FROM_VAF {
    tag "$meta.id"
    label 'process_single'
    label 'memory_medium'

    container 'docker.io/ferranmuinos/test_mutated_genomes'

    input:
    tuple val(meta), path (mutated_genomes_results)
    path (clinical_features)

    output:
    tuple val(meta), path("*.tsv") , emit: mutated_cells_sample
    path  "versions.yml" , topic: versions


    script:
    """
    mutgenomes_summary_tables.py \\
        --metadata-file ${clinical_features} ;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    """
    touch ${prefix}.sample.mutated_genomes_from_vaf.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}