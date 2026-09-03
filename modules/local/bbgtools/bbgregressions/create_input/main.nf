process CREATE_INPUT {

    tag "regressions"

    label "bbgregressions"

    input:
    path config
    path data

    output:
    path ("input/*")        , emit:  inputs
    path "versions.yml"     , topic: versions

    script:
    """
    bbgregressions create_input -config ${config};

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
