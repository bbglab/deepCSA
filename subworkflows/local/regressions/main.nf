include { RUNREGRESSIONS    as RUNREGRESSIONS             } from '../../../modules/local/runregressions/main'

workflow REGRESSIONS{

    take:
    metric_name
    metric_data
    metric_params

    main:
    ch_versions = Channel.empty()

    RUNREGRESSIONS(metric_name,
    metric_data, metric_params,
    params.responses_subset_regressions, params.responses_excluded_regressions,
    params.samples_subset_regressions,
    params.predictors_file_regressions, params.predictors_plot_config_regressions,
    params.random_effects_vars_regressions, params.multipletesting_join_regressions,
    params.multivariate_rules_regressions, params.response_subplots,
    params.total_plot, params.response_and_total_subplots, params.make2,
    params.correct_pvals, params.sign_threshold)
    ch_versions = ch_versions.mix(RUNREGRESSIONS.out.versions)

    emit:
    res_pdf = RUNREGRESSIONS.out.res_pdf
    res_tables = RUNREGRESSIONS.out.res_tables
    versions = ch_versions                // channel: [ versions.yml ]
}
