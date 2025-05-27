
process NORMALIZE_SIGNATURES {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/ferriolcalvet/hdp_stefano:0.1.0'


    input:
    tuple val(meta), path(output_dir)

    output:
    tuple val(meta), path("normalized_output_dir"), emit: normalized_results
    path "versions.yml"                           , topic: versions


    script:
    """
    Rscript run_HDP_sigNormalising.R \\
        $output_dir \\
        ${params.norm_file}
    cp -r $output_dir normalized_output_dir

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


