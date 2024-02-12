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
    // TODO see if we provide directly the panel annotated as omega requires
    // see regions_sites.annotation_summary.tsv in some omega test directory

    output:
    tuple val(meta), path("mutability_per_sample_gene_context.tsv"), path("mutations_per_sample_gene_impact_context.tsv") , emit: mutabs_n_mutations_tsv
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    omega preprocessing --preprocessing-mode compute_mutabilities \\
                        --depths-file ${depths} \\
                        --mutations-file ${mutations} \\
                        --input-vep-postprocessed-file ${annotated_panel} \\
                        --mutabilities-table mutability_per_sample_gene_context.tsv \\
                        --table-observed-muts mutations_per_sample_gene_impact_context.tsv
    #--mutational-profile ${mutation_profile}
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

