process CREATE_INPUT {

    tag "regressions"
    label 'process_single'

    container "docker.io/bbglab/bbgregressions:0.1.0"

    input:
    path config
    path data

    output:
    path ("inputs/*")       , emit:  inputs
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
