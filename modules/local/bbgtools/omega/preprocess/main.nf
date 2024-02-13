process OMEGA_PREPROCESS {
    tag "$meta.id"
    label 'process_high'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'docker.io/ferriolcalvet/omega:latest'

    input:
    tuple val(meta), path(mutations), path(depths), path(mutation_profile)
    tuple val(meta2), path (annotated_panel)

    output:
    tuple val(meta), path("mutability_per_sample_gene_context.*.tsv"), path("mutations_per_sample_gene_impact_context.*.tsv") , emit: mutabs_n_mutations_tsv
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO revise this fix
    def sample_name = prefix.tokenize('.')[0]
    """
    omega preprocessing --preprocessing-mode compute_mutabilities \\
                        --depths-file ${depths} \\
                        --mutations-file ${mutations} \\
                        --input-vep-postprocessed-file ${annotated_panel} \\
                        --mutabilities-table mutability_per_sample_gene_context.${prefix}.tsv \\
                        --table-observed-muts mutations_per_sample_gene_impact_context.${prefix}.tsv \\
                        --mutational-profile ${mutation_profile} \\
                        --single-sample ${sample_name}
    # $args -c $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omega: 1.0
    END_VERSIONS
    """
}

