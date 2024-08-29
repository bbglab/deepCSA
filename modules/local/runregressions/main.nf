process RUNREGRESSIONS {

    tag "regressions"
    label 'process_single'

    container 'docker.io/rblancomi/statsmodels:test'

    input:
    path(all_mutrates_file)
    path(oncodrivefml_regressions_files)
    path(omega_regressions_files)
    val(mutrate_regressions)
    val(omega_regressions)
    val(oncodrivefml_regressions)
    val(responses_subset_regressions)
    val(samples_subset_regressions)
    val(predictors_file_regressions)
    val(predictors_plot_config_regressions)
    val(random_effects_vars_regressions)
    val(multipletesting_join_regressions)
    val(multivariate_rules_regressions)
    val(response_subplots)
    val(total_plot)
    val(response_and_total_subplots)
    val(make2)
    val(correct_pvals)
    val(sign_threshold)

    output:
    path("regressions.pdf")            , emit: regressions
    path  "versions.yml"               , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "regressions"
    def predictors_plot_config = "\'${groovy.json.JsonOutput.toJson(predictors_plot_config_regressions)}\'"
    def multipletesting_join = "\'${groovy.json.JsonOutput.toJson(multipletesting_join_regressions)}\'"
    def multivariate_rules = "\'${groovy.json.JsonOutput.toJson(multivariate_rules_regressions)}\'"
    def oncodrivefml_regressions_files_str = oncodrivefml_regressions_files.join(',')
    def omega_regressions_files_str = omega_regressions_files.join(',')

    """
    regression_analysis.py \\
                -pdf regressions.pdf \\
                --mutrate_data ${all_mutrates_file} \\
                --oncodrivefml_data ${oncodrivefml_regressions_files_str} \\
                --omega_data ${omega_regressions_files_str} \\
                --mutrate_params ${mutrate_regressions} \\
                --omega_params ${omega_regressions} \\
                --oncodrivefml_params ${oncodrivefml_regressions} \\
                --responses_subset ${responses_subset_regressions} \\
                --samples_subset ${samples_subset_regressions} \\
                --predictors_file ${predictors_file_regressions} \\
                --predictors_plot_config ${predictors_plot_config} \\
                --random_effects_vars ${random_effects_vars_regressions} \\
                --multipletesting_join ${multipletesting_join} \\
                --multivariate_rules ${multivariate_rules} \\
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
    def prefix = task.ext.prefix ?: "regressions"
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
