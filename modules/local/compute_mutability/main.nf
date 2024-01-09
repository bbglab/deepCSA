process COMPUTE_MUTABILITY {

    tag "$meta.id"
    label 'process_single'
    label 'process_high_memory'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'docker.io/ferriolcalvet/bgreference'

    input:
    tuple val(meta), path(matrix), path(mut_profile)
    tuple val(meta2), path(depths)

    output:
    tuple val(meta), path("*.mutability_per_site.tsv") , emit: mutability
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mutprof_3compute_mutability.py \\
                    --sample_name ${prefix} \\
                    --mutation_matrix ${matrix} \\
                    --depths ${depths} \\
                    --profile ${mut_profile} \\
                    --bedfile ${params.bedf} \\
                    --out_mutability ${prefix}.mutability_per_site.tsv \\
                    ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """


    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.profile.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
