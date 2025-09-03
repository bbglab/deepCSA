////
// THIS SCRIPT IS CURRENTLY NOT IN USE
////


process SIGPROFILEREXTRACTOR {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/ferriolcalvet/sigprofilerassignment'

    input:
    tuple val(meta), path(matrix)

    output:
    tuple val(meta), path("**.pdf")                                         , emit: plots
    tuple val(meta), path("**.txt")                                         , emit: stats
    tuple val(meta), path("**Decomposed_MutationType_Probabilities.*.txt")  , emit: mutation_probs
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def assembly = task.ext.assembly ?: "GRCh38"
    // python -c "from SigProfilerAssignment import Analyzer as Analyze; Analyze.cosmic_fit('${matrix}', 'output', input_type='matrix', context_type='96',
    //                    collapse_to_SBS96=True, cosmic_version=3.4, exome=False,
    //                    genome_build="GRCh38", signature_database='${reference_signatures}',
    //                    exclude_signature_subgroups=None, export_probabilities=False,
    //                    export_probabilities_per_mutation=False, make_plots=True,
    //                    sample_reconstruction_plots=False, verbose=False)"
    // mv output_${prefix}/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities.txt output_${prefix}/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities.${prefix}.txt;
    """
    signatures_sigprofilerextractor.py  matrix . ${matrix} \\
                    --ref_genome ${assembly} \\
                    ${args} \\
                    --cpu ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        SigProfilerAssignment : 0.1.1
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
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
