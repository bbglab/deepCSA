process OMEGA_ESTIMATOR {
    tag "$meta.id"
    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'


    container 'docker.io/ferriolcalvet/omega:latest'

    input:
    tuple val(meta) , path(mutabilities_table), path(mutations_table), path(depths)
    tuple val(meta2), path(annotated_panel)
    path (genes_json)

    output:
    tuple val(meta), path("output_*.tsv"), emit: results
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def option = task.ext.option ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    mkdir groups;

    mv ${genes_json} groups/group_genes.json

    cat > groups/group_impacts.json << EOF
    {
        "missense": ["missense"],
        "nonsense": ["nonsense"],
        "essential_splice": ["essential_splice"],
        "splice_region_variant": ["splice_region_variant"],
        "truncating": ["nonsense", "essential_splice"],
        "essential_splice_plus": ["essential_splice", "splice_region_variant"],
        "truncating_plus": ["nonsense", "essential_splice", "splice_region_variant"],
        "nonsynonymous_splice": ["missense", "nonsense", "essential_splice"]
    }
    EOF

    cat > groups/group_samples.json << EOF
    {
        "${meta.id}" : ["${meta.id}"]
    }
    EOF

    omega estimator --mutability-file ${mutabilities_table} \\
                    --observed-mutations-file ${mutations_table} \\
                    --depths-file ${depths} \\
                    --vep-annotation-file ${annotated_panel} \\
                    --grouping-folder ./groups/ \\
                    --output-fn output_${option}.${prefix}.tsv \\
                    --option ${option} \\
                    --cores 4
                    # --cores ${task.cpus}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def option = task.ext.option ?: "bayes"
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch output_${option}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """
}

