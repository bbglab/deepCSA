process OMEGA_ESTIMATOR {
    tag "$meta.id"
    label 'cpu_medium'


    container 'docker.io/bbglab/omega:0.2.1'

    input:
    tuple val(meta) , path(mutabilities_table), path(mutations_table), path(depths)
    tuple val(meta2), path(annotated_panel)
    path (genes_json)
    path (impacts_json)

    output:
    tuple val(meta), path("output_*.tsv"), emit: results
    path "versions.yml"                  , topic: versions


    script:
    def option = task.ext.option ?: ""
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    mkdir groups;

    mv ${genes_json} groups/group_genes.json
    mv ${impacts_json} groups/group_impacts.json
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
                    --cores ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: \$(omega --version | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    def option = task.ext.option ?: ""
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch output_${option}.${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: \$(omega --version | cut -d ' ' -f 3)
    END_VERSIONS
    """
}

