process ANNOTATE_OMEGA_QC {

    tag "all_samples"
    label 'process_low'

    label 'deepcsa_core'

    input:
    path (all_omegas)
    path (compiled_flagged_cases)

    output:
    path("*flagged_annotated.tsv")                      , emit: all_omegas_annotated
    path("debug.*flagged*.tsv")         , optional: true, emit: flagged_cases
    path("debug.syn_flagged_gene.tsv")  , optional: true, emit: flagged_synonymous_cases
    path("*.tsv")                       , optional: true, emit: files
    path("*.png")                       , optional: true, emit: plots
    path "versions.yml"                 , topic: versions

    script:
    """
    ls ${compiled_flagged_cases} > compiled_flagging_cases.txt;

    annotate_omega_failing.py \\
                    --omegas-file ${all_omegas} \\
                    --compiled-flagged-files compiled_flagging_cases.txt \\
                    --output omega.flagged_annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "all_samples"
    """
    touch ${prefix}.pdf
    touch omega.flagged_annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
