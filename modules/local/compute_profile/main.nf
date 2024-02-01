process COMPUTE_PROFILE {

    tag "$meta.id"
    label 'process_low'

    // // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(matrix), path(trinucleotide)

    output:
    tuple val(meta), path("*.profile.tsv")                             , emit: profile
    tuple val(meta), path("*.pdf")                    , optional:true  , emit: plots
    tuple val(meta), path("*.matrix.WGS")             , optional:true  , emit: wgs
    tuple val(meta), path("*.matrix.WGS.sigprofiler") , optional:true  , emit: wgs_sigprofiler
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filters = task.ext.filters ?: ""

    """
    mut_profile.py profile \\
                    --sample_name ${prefix} \\
                    --mutation_matrix ${matrix} \\
                    --trinucleotide_counts ${trinucleotide} \\
                    --out_profile ${prefix}.profile.tsv \\
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
