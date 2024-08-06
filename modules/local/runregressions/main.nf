process RUNREGRESSIONS {

    tag "regressions"
    label 'process_single'

    container 'docker.io/rblancomi/statsmodels:test'

    input:
    path(regressions_config)

    output:
    path  "regressions.pdf"             , emit: regressions
    path  "versions.yml"                , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "regressions"

    """
    regression_analysis.py \\
                -config ${regressions_config} \\
                -pdf regressions.pdf \\
                --response_subplots \\
                --no-total_plot \\
                --no-response_and_total_subplots \\
                --make2 \\
                --correct_pvals \\
                --sign_threshold 0.2;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "regression"
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
