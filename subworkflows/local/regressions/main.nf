include { RUNREGRESSIONS    as RUNREGRESSIONS             } from '../../../modules/local/runregressions/main'

workflow REGRESSIONS{
    take:
    all_mutrates
    oncodrivefml_regressions_files
    omega_regressions_files
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

    main:
    ch_versions = Channel.empty()

    RUNREGRESSIONS(all_mutrates, oncodrivefml_regressions_files,
    omega_regressions_files, mutrate_regressions,
    omega_regressions, oncodrivefml_regressions,
    responses_subset_regressions, samples_subset_regressions,
    predictors_file_regressions, predictors_plot_config_regressions,
    random_effects_vars_regressions, multipletesting_join_regressions,
    multivariate_rules_regressions, response_subplots,
    total_plot, response_and_total_subplots, make2,
    correct_pvals, sign_threshold)
    ch_versions = ch_versions.mix(RUNREGRESSIONS.out.versions)

    emit:
    regressions = RUNREGRESSIONS.out.regressions
    versions = ch_versions                // channel: [ versions.yml ]
}
