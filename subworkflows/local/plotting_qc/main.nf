
include { PLOT_MUTDENSITY_QC        as PLOTMUTDENSITYQC     } from '../../../modules/local/plot/qc/mutation_densities/main'
include { PLOT_METRICS_VS_DEPTH_QC  as PLOTMETRICSVSDEPTHQC } from '../../../modules/local/plot/qc/metrics_vs_depth/main'
include { ANNOTATE_OMEGA_QC         as APPLYOMEGAQC         } from '../../../modules/local/plot/qc/annotate_omega/main'
include { PLOT_MUTATION_SPECIFIC    as PLOTMUTATIONSPECIFIC } from '../../../modules/local/plot/qc/mutation_specific/main'
include { PLOT_OMEGA_VS_DNDSCV      as PLOTOMEGAVSDNDSCV    } from '../../../modules/local/plot/qc/omega_vs_dndscv/main'


workflow PLOTTING_QC {

    take:
    all_mutations
    all_mutdensities
    all_adjusted_mutdensities
    all_omegas_globalloc
    average_depth_gene_sample
    all_omegas
    panel
    groups_definition
    group_name
    dndscv_cv 
    // all_samples_depth
    // all_groups
    // full_panel_rich
    // seqinfo_df
    // domain_df
    // exons_depths_df


    main:

    // Channel.of([ [ id: "all_samples" ] ])
    // .join( all_mutations )
    // .set{ mutations }

    // plotting only for the entire cohort group
    // channel.of([ [ id: "all_samples" ] ])
    // .join( positive_selection_results_ready )
    // .set{ all_samples_results }

    PLOTMUTATIONSPECIFIC(all_mutations)

    PLOTMUTDENSITYQC(all_mutdensities, panel, groups_definition, group_name)

    PLOTMETRICSVSDEPTHQC(
        all_mutdensities,
        average_depth_gene_sample.map { it -> it[1] },
        groups_definition,
        group_name,
        all_adjusted_mutdensities,
        all_omegas_globalloc
    )


    APPLYOMEGAQC(all_omegas, PLOTMUTDENSITYQC.out.compiled_flagged.collect())
    
    // Only run omega vs dndscv QC plot when dndscv results are available
    omega_vs_dndscv_plots = channel.empty()
    if (params.dnds) {
        PLOTOMEGAVSDNDSCV(all_omegas, dndscv_cv, groups_definition, group_name, APPLYOMEGAQC.out.flagged_cases)
        omega_vs_dndscv_plots = PLOTOMEGAVSDNDSCV.out.plots
    }
    

    emit:
    mutdensity_plots        = PLOTMUTDENSITYQC.out.plots
    metrics_vs_depth_plots  = PLOTMETRICSVSDEPTHQC.out.plots
    metrics_vs_depth_tables = PLOTMETRICSVSDEPTHQC.out.tables
    flagged_omegas          = APPLYOMEGAQC.out.all_omegas_annotated
    omega_vs_dndscv_qc      = omega_vs_dndscv_plots

}
