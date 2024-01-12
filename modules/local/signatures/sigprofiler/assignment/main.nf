process SIGPROFILERASSIGNMENT {
    tag "$meta.id"
    label 'process_medium'


    container 'docker.io/ferriolcalvet/sigprofilerassignment'

    input:
    tuple val(meta), path(matrix)
    path(reference_signatures)

    output:
    tuple val(meta), path("**.pdf"), emit: plots
    tuple val(meta), path("**.txt"), emit: stats
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // python -c "from SigProfilerAssignment import Analyzer as Analyze; Analyze.cosmic_fit('${matrix}', 'output', input_type='matrix', context_type='96',
    //                    collapse_to_SBS96=True, cosmic_version=3.4, exome=False,
    //                    genome_build="GRCh38", signature_database='${reference_signatures}',
    //                    exclude_signature_subgroups=None, export_probabilities=False,
    //                    export_probabilities_per_mutation=False, make_plots=True,
    //                    sample_reconstruction_plots=False, verbose=False)"
    """
    #python -c "from SigProfilerAssignment import Analyzer as Analyze; Analyze.cosmic_fit('${matrix}', 'output', input_type='matrix', context_type='96', signature_database='${reference_signatures}', genome_build='GRCh38', sample_reconstruction_plots= 'pdf', exclude_signature_subgroups=['Chemotherapy_signatures','Immunosuppressants_signatures'])"
    python -c "from SigProfilerAssignment import Analyzer as Analyze; Analyze.cosmic_fit('${matrix}', 'output', input_type='matrix', context_type='96', genome_build='GRCh38', sample_reconstruction_plots= 'pdf', exclude_signature_subgroups=['Chemotherapy_signatures','Immunosuppressants_signatures'])"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
