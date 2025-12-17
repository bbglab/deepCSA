process REGRESSIONS {

    tag "regressions"
    label 'process_single'

    container "docker.io/bbglab/bbgregressions:0.1.0"

    input:
    path config

    output:
    path ("regressions/*")       , emit: models
    path "versions.yml"          , topic: versions

    script:
    """
    bbgregressions regressions -config ${config};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
