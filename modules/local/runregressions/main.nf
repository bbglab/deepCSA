process RUNREGRESSIONS {

    tag "regressions"
    label 'process_single'

    container 'docker.io/rblancomi/statsmodels:test'

    input:
    path(features_table)

    output:
    path("regressions.pdf")             , emit: regressions
    path  "versions.yml"                , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "regressions"

    """
    echo "Starting regression analysis..."
    echo "Current directory: \$(pwd)"
    echo "Listing files in current directory:"
    ls -l
    ls ../../../bin/regression_analysis.py

    regression_analysis.py \\
                -config /workspace/datasets/transfer/ferriol_deepcsa/regressions/config.ini \\
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
