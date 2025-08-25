process COMPARE_SIGNATURES {
    tag "${meta.id}"
    label 'process_medium'

    container 'docker.io/ferriolcalvet/hdp_stefano:0.1.0'

    input:
    tuple val(meta), path(output_dir)
    path reference_signatures

    output:
    tuple val(meta), path("*.compared_output_dir/**"), emit: compared_results
    path "versions.yml", topic: versions

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    Rscript /app/HDP_sigExtraction/R/run_HDP_comparing.R \\
        output_dir/ \\
        ${args} \\
        ${reference_signatures} \\
        "NA" \\
        "FALSE"
    cp -r output_dir ${prefix}.compared_output_dir

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R       : \$(Rscript --version | sed -e 's/.*version //g')
        Rscript : \$(Rscript --version | sed -e 's/.*version //g')
        HDP     : original
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R       : \$(Rscript --version | sed -e 's/.*version //g')
        Rscript : \$(Rscript --version | sed -e 's/.*version //g')
        HDP     : original
    END_VERSIONS
    """
}
