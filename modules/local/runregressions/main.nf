process RUNREGRESSIONS {

    tag "regressions"
    label 'process_single'

    container 'docker.io/rblancomi/statsmodels:test'

    input:
    all_mutrates
    path oncodrivefml_regressions_files_dir from oncodrivefml_regressions_files.collect()
    path omega_regressions_files_dir from omega_regressions_files.collect()
    mutrate_regressions
    omega_regressions
    oncodrivefml_regressions
    responses_subset_regressions
    samples_subset_regressions
    predictors_file_regressions
    predictors_plot_config_regressions
    random_effects_vars_regressions
    multipletesting_join_regressions
    multivariate_rules_regressions
    response_subplots
    total_plot
    response_and_total_subplots
    make2
    correct_pvals
    sign_threshold

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
                -pdf regressions.pdf \\
                --mutrate_file ${all_mutrates} \\
                --oncodrivefml_regressions_dir ${oncodrivefml_regressions_files_dir} \\
                --omega_regressions_dir ${omega_regressions_files_dir} \\
                --mutrate_params ${mutrate_regressions} \\
                --omega_params ${omega_regressions} \\
                --oncodrivefml_params ${oncodrivefml_regressions} \\
                --responses_subset ${responses_subset_regressions} \\
                --samples_subset ${samples_subset_regressions} \\
                --predictors_file ${predictors_file_regressions} \\
                --predictors_plot_config ${predictors_plot_config_regressions} \\
                --random_effects_vars ${random_effects_vars_regressions} \\
                --multipletesting_join ${multipletesting_join_regressions} \\
                --multivariate_rules ${multivariate_rules_regressions} \\
                --response_subplots ${response_subplots} \\
                --total_plot ${total_plot} \\
                --response_and_total_subplots ${response_and_total_subplots} \\
                --make2 ${make2} \\
                --correct_pvals ${correct_pvals} \\
                --sign_threshold ${sign_threshold};

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
