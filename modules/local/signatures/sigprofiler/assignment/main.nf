process SIGPROFILERASSIGNMENT {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/ferriolcalvet/sigprofilerassignment'

    input:
    tuple val(meta), val(type), path(matrix)
    path(reference_signatures)

    output:
    tuple val(meta), path("**.pdf")                                         , emit: plots
    tuple val(meta), path("**.txt")                                         , emit: stats
    tuple val(meta), path("**Decomposed_MutationType_Probabilities.*.txt")  , emit: mutation_probs
    path "versions.yml"                                                     , topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${type}"
    def assembly = task.ext.assembly ?: "GRCh38"
    // python -c "from SigProfilerAssignment import Analyzer as Analyze; Analyze.cosmic_fit('${matrix}', 'output', input_type='matrix', context_type='96',
    //                    collapse_to_SBS96=True, cosmic_version=3.4, exome=False,
    //                    genome_build="GRCh38", signature_database='${reference_signatures}',
    //                    exclude_signature_subgroups=None, export_probabilities=False,
    //                    export_probabilities_per_mutation=False, make_plots=True,
    //                    sample_reconstruction_plots=False, verbose=False)"
    """
    #python -c "from SigProfilerAssignment import Analyzer as Analyze; Analyze.cosmic_fit('${matrix}', 'output_${prefix}', input_type='matrix', context_type='96', signature_database='${reference_signatures}', genome_build='${assembly}', sample_reconstruction_plots= 'pdf', exclude_signature_subgroups= ${params.exclude_subgroups})"
    python -c "from SigProfilerAssignment import Analyzer as Analyze; Analyze.cosmic_fit('${matrix}', 'output_${prefix}', input_type='matrix', context_type='96', genome_build='${assembly}', signature_database='${reference_signatures}', exclude_signature_subgroups=${params.exclude_subgroups})"

    mv output_${prefix}/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities.txt output_${prefix}/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities.${prefix}.txt;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        SigProfilerAssignment : 0.1.1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        SigProfilerAssignment : 0.1.1
    END_VERSIONS
    """
}
