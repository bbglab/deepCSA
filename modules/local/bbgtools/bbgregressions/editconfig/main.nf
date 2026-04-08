process EDITCONFIG {

    tag "regressions"
    label 'process_low'

    container "docker.io/rblancomi/bbgregressions:dev"

    input:
    path config
    val mode
    val metric
    path omega_res
    path groups

    output:
    path ("config.yaml")        , emit:  config
    path "versions.yml"                 , topic: versions

    script:
    """
    regressions_editconfig.py \\
                        --config_file ${config} \\
                        --mode ${mode} \\
                        --metric ${metric} \\
                        --omega_res_file ${omega_res} \\
                        --groups_file ${groups} \\

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
