include { RUNREGRESSIONS    as RUNREGRESSIONS             } from '../../../modules/local/runregressions/main'

workflow REGRESSIONS{

    take:
    metric_name
    metric_data
    metric_params

    main:

    predictors_file = file(params.predictors_file_regressions)


    RUNREGRESSIONS(metric_name,
                    metric_data, metric_params,
                    params.responses_subset_regressions, params.responses_excluded_regressions,
                    params.samples_subset_regressions,
                    predictors_file, params.predictors_plot_config_regressions,
                    params.random_effects_vars_regressions, params.multipletesting_join_regressions,
                    params.multivariate_rules_regressions, params.response_subplots
                    )


    emit:
    res_pdf = RUNREGRESSIONS.out.res_pdf
    res_tables = RUNREGRESSIONS.out.res_tables

}
