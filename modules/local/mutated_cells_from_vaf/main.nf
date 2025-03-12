process MUTATED_CELLS_FROM_VAF {
    tag "$meta.id"
    label 'process_single'
    label 'memory_medium'

    container 'docker.io/ferranmuinos/test_mutated_genomes'
    // conda activate probabilistic

    input:
    tuple val(meta), path (mutated_genomes_results)
    path (clinical_features)

    output:
    tuple val(meta), path("*.tsv") , emit: mutated_cells_sample
    path  "versions.yml" , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // def separator = task.ext.separator ?: "comma"
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
    // def panel_version = task.ext.panel_version ?: "${meta2.id}"
    """
    touch ${prefix}.sample.mutated_genomes_from_vaf.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}

// this could be used to define the way sex is encodedand which is the sample id column
//     def features = task.ext.features ?: ""
//     """
//     cat > features_table_information.json << EOF
//     {
//         ${features}
//     }
//     EOF
