process COMPUTE_MUTABILITY {

    tag "$meta.id"
    label 'process_low_fixed_cpus'
    label 'process_high_memory'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta) , path(matrix), path(mut_profile), path(depths)
    tuple val(meta2), path(bedfile)

    output:
    // TODO revise this to see which one is outputed and why
    tuple val(meta), path("*.mutability_per_site.tsv")                           , emit: mutability_not_adjusted
    tuple val(meta), path("*.mutability_per_site.tsv.adjusted")                  , emit: mutability
    path "versions.yml"                                                          , topic: versions

    // tuple val(meta), path("*.mutability_per_site.tsv")                           , emit: mutability
    // tuple val(meta), path("*.mutability_per_site.tsv.adjusted") , optional:true  , emit: mutability_adjusted
    // path "versions.yml"                                                          , topic: versions

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
                    --bedfile ${bedfile} \\
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
