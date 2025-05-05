process OMEGA_MUTABILITIES {
    tag "$meta.id"
    label 'cpu_single_fixed'
    label 'time_low'
    label 'process_high_memory'


    container 'docker.io/ferriolcalvet/omega:20250121'

    input:
    tuple val(meta) , path(mutabilities_table), path(mutations_table), path(depths)
    tuple val(meta2), path(annotated_panel)
    path (genes_json)

    output:
    tuple val(meta), path("mutabilities_per_site.*.tsv.gz") , emit: mutabilities
    path "versions.yml"                                     , topic: versions


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    mkdir groups;

    mv ${genes_json} groups/group_genes.json

    cat > groups/group_samples.json << EOF
    {
        "${meta.id}" : ["${meta.id}"]
    }
    EOF

    omega mutabilities --mutability-file ${mutabilities_table} \\
                        --depths-file ${depths} \\
                        --vep-annotation-file ${annotated_panel} \\
                        --grouping-folder ./groups/ \\
                        --output-fn mutabilities_per_site.${prefix}.tsv.gz \\
                        --cores 4
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch mutabilities_per_site.${prefix}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """
}

