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
    def name = "${meta.id}.${type}"
    def assembly = task.ext.assembly ?: "GRCh38"
    """
    python -c "from SigProfilerAssignment import Analyzer as Analyze; Analyze.cosmic_fit('${matrix}', 'output_${name}', input_type='matrix', context_type='96', genome_build='${assembly}', signature_database='${reference_signatures}', exclude_signature_subgroups=${params.exclude_subgroups})"

    mv output_${name}/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities.txt output_${name}/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities.${name}.txt;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        SigProfilerAssignment : 0.1.1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        SigProfilerAssignment : 0.1.1
    END_VERSIONS
    """
}
