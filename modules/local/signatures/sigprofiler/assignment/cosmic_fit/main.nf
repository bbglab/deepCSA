process SIGPROFILERASSIGNMENT_COSMIC_FIT {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/ferriolcalvet/sigprofiler_assignment:1.1.3'

    input:
    tuple val(meta), val(type), path(matrix)
    path(reference_signatures)

    output:
    tuple val(meta), path("**.pdf")                                         , emit: plots
    tuple val(meta), path("**.txt")                                         , emit: stats
    tuple val(meta), path("**Decomposed_MutationType_Probabilities.*.txt")  , emit: mutation_probs
    path "versions.yml"                                                     , topic: versions


    script:
    def assembly = task.ext.genome_assembly ?: "GRCh38"
    def context_type = task.ext.context_type ?: "96" // "96", "288", "1536", "DINUC", and "ID"
    def name = "${meta.id}.${type}.${context_type}"

    // FIXME: the definition of subgroups to exclude seems not to work in the new CLI SigProfilerAssignment
    // def exclude_signature_subgroups = params.exclude_subgroups ? "--exclude_signature_subgroups \"${params.exclude_subgroups}\"" : ""
    // The default input format is "matrix".
    """
    mkdir -p spa_volume
    export SIGPROFILERMATRIXGENERATOR_VOLUME='./spa_volume'
    export SIGPROFILERPLOTTING_VOLUME='./spa_volume'
    export SIGPROFILERASSIGNMENT_VOLUME='./spa_volume'

    SigProfilerAssignment cosmic_fit \\
            ${matrix} \\    
            output_${name} \\
            --signature_database ${reference_signatures} \\
            --genome_build ${assembly} \\
            --cpu ${task.cpus} \\
            --context_type ${context_type} \\
            --volume spa_volume

    cp output_${name}/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities.txt output_${name}/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities.${name}.txt;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
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
        SigProfilerAssignment : 1.1.3
    END_VERSIONS
    """
}
