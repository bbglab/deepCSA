include { RUNREGRESSIONS    as RUNREGRESSIONSCOHORT             } from '../../../modules/local/runregressions/main'
include { RUNREGRESSIONS    as RUNREGRESSIONSIGNORE             } from '../../../modules/local/runregressions/main'
include { RUNREGRESSIONS    as RUNREGRESSIONSMEAN               } from '../../../modules/local/runregressions/main'


workflow REGRESSIONS{

    take:
    metric_name
    metric_data
    metric_params

    main:

    predictors_file = file(params.predictors_file_regressions)


    RUNREGRESSIONSCOHORT(metric_name,
                    metric_data, metric_params,
                    params.responses_subset_regressions, params.responses_excluded_regressions,
                    params.samples_subset_regressions,
                    predictors_file, params.predictors_plot_config_regressions,
                    params.multipletesting_join_regressions,
                    params.multivariate_rules_regressions
                    )


    RUNREGRESSIONSIGNORE(metric_name,
                    metric_data, metric_params,
                    params.responses_subset_regressions, params.responses_excluded_regressions,
                    params.samples_subset_regressions,
                    predictors_file, params.predictors_plot_config_regressions,
                    params.multipletesting_join_regressions,
                    params.multivariate_rules_regressions
                    )


    RUNREGRESSIONSMEAN(metric_name,
                    metric_data, metric_params,
                    params.responses_subset_regressions, params.responses_excluded_regressions,
                    params.samples_subset_regressions,
                    predictors_file, params.predictors_plot_config_regressions,
                    params.multipletesting_join_regressions,
                    params.multivariate_rules_regressions
                    )


    emit:
    res_pdf = RUNREGRESSIONSCOHORT.out.res_pdf
    res_tables = RUNREGRESSIONSCOHORT.out.res_tables

}
