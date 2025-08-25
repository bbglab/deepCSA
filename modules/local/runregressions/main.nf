process RUNREGRESSIONS {

    tag "regressions"
    label 'process_single'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    val metric_name
    path metric_data
    val metric_params

    val responses_subset_regressions
    val responses_excluded_regressions
    val samples_subset_regressions
    path predictors_file_regressions
    val predictors_plot_config_regressions
    val multipletesting_join_regressions
    val multivariate_rules_regressions

    output:
    path ("inputs/*"), emit: regression_inputs
    path (metric_name), emit: res_tables
    path ("*.pdf"), emit: res_pdf
    path "versions.yml", topic: versions

    script:
    def args = task.ext.args ?: ""
    def responses_excluded = "\'${groovy.json.JsonOutput.toJson(responses_excluded_regressions)}\'"
    def predictors_plot_config = "\'${groovy.json.JsonOutput.toJson(predictors_plot_config_regressions)}\'"
    def multipletesting_join = "\'${groovy.json.JsonOutput.toJson(multipletesting_join_regressions)}\'"
    def multivariate_rules = "\'${groovy.json.JsonOutput.toJson(multivariate_rules_regressions)}\'"
    def metric_data_str = metric_data.join(',')
    """
    regression_analysis.py \\
                -pdf ${metric_name}.pdf \\
                --metric ${metric_name} \\
                --data ${metric_data_str} \\
                --predictors_file ${predictors_file_regressions} \\
                --metric_params ${metric_params} \\
                --responses_subset ${responses_subset_regressions} \\
                --responses_excluded ${responses_excluded} \\
                --samples_subset ${samples_subset_regressions} \\
                --predictors_plot_config ${predictors_plot_config} \\
                --multipletesting_join ${multipletesting_join} \\
                --multivariate_rules ${multivariate_rules} \\
                --save_tables_dir ${metric_name} \\
                ${args};

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
