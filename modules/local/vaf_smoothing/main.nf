process VAF_SMOOTHING {
    tag "groups"
    label 'process_low'

    label 'deepcsa_core'

    input:
    tuple val(meta) , path(all_mutations) 
    path (all_mutdensities) 
    tuple val(meta2), path(average_depth_sample) 
    
    output:
    path("*.tsv.gz")                        , emit: smoothed_vaf_tables
    path("*.pdf")           , optional: true, emit: plots
    path  "versions.yml"                    , topic: versions

    script:
    """
    vaf_smoothing.py \\
            --mutations ${all_mutations} \\
            --mutdensities ${all_mutdensities} \\
            --depth-sample ${average_depth_sample}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch vaf_smoothing.pdf
    touch smoothed_VAFs.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}