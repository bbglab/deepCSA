process SIGPROFILERASSIGNMENT_DECOMPOSE_FIT {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/ferriolcalvet/sigprofiler_assignment:1.1.3'

    input:
    tuple val(meta), val(type), path(matrix)
    path (extracted_signatures)
    path (reference_signatures)

    output:
    tuple val(meta), path("**.pdf")                                         , emit: plots
    tuple val(meta), path("**.txt")                                         , emit: stats
    tuple val(meta), path("**Decomposed_MutationType_Probabilities.*.txt")  , emit: mutation_probs
    path "versions.yml"                                                     , topic: versions


    script:
    def name = "${meta.id}.${type}"
    def assembly = task.ext.assembly ?: "GRCh38"
    
    // FIXME: the definition of subgroups to exclude seems not to work in the new CLI SigProfilerAssignment
    // def exclude_signature_subgroups = params.exclude_subgroups ? "--exclude_signature_subgroups \"${params.exclude_subgroups}\"" : ""
    """
    mkdir -p spa_volume
    export SIGPROFILERMATRIXGENERATOR_VOLUME='./spa_volume'
    export SIGPROFILERPLOTTING_VOLUME='./spa_volume'
    export SIGPROFILERASSIGNMENT_VOLUME='./spa_volume'

    SigProfilerAssignment decompose_fit \\
            ${matrix} \\
            signature_decomposition_${name} \\
            --signatures ${extracted_signatures} \\
            --signature_database ${reference_signatures} \\
            --genome_build ${assembly} \\
            --cpu ${task.cpus} \\
            --volume spa_volume

    #mv signature_decomposition_${name}/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities.txt signature_decomposition_${name}/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities.${name}.txt;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        SigProfilerAssignment : 1.1.3
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
        SigProfilerAssignment : 1.1.3
    END_VERSIONS
    """
}
