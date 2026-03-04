process SIGPROFILER_MATRIXGENERATOR {
    tag "samples"
    label 'process_single'

    container 'docker.io/ferriolcalvet/sigprofilermatrixgenerator:1.3.5'

    input:
    path (vcf)

    output:
    path("**plots/*"), optional : true, emit: output_plots
    path("**ID/*")   , optional : true, emit: matrices_ID
    path("**DBS/*")  , optional : true, emit: matrices_DBS
    path("**SBS/*")  , optional : true, emit: matrices_SBS
    path("**TSB/*")  , optional : true, emit: transcription_bias
    path "versions.yml"                                    , topic: versions


    script:
    def prefix = task.ext.prefix ?: 'samples'
    def args = task.ext.args ?: ""
    def genome = task.ext.genome_assembly ?: "GRCh38"
    """
    mkdir input_mutations
    cp *.vcf input_mutations/.

    SigProfilerMatrixGenerator matrix_generator \\
                ${prefix} \\
                ${genome} \\
                input_mutations/ \\
                ${args}
            
    mv input_mutations/output/ ${prefix}_output/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SigProfilerMatrixGenerator: 1.3.5
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SigProfilerMatrixGenerator: 1.3.5
    END_VERSIONS
    """
}