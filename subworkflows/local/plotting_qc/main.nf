
include { PLOT_MUTDENSITY_QC    as PLOTMUTDENSITYQC     } from '../../../modules/local/plot/qc/mutation_densities/main'
include { ANNOTATE_OMEGA_QC     as APPLYOMEGAQC         } from '../../../modules/local/plot/qc/annotate_omega/main'


workflow PLOTTING_QC {

    take:
    // positive_selection_results_ready
    all_mutdensities
    // all_samples_depth
    // all_groups
    all_omegas
    panel
    groups_definition
    group_name
    // full_panel_rich
    // seqinfo_df
    // domain_df
    // exons_depths_df


    main:

    // pdb_tool_df   = params.annotations3d
    //                         ? channel.fromPath( "${params.annotations3d}/pdb_tool_df.tsv", checkIfExists: true).first()
                            // : channel.empty()


    // plotting only for the entire cohort group
    // channel.of([ [ id: "all_samples" ] ])
    // .join( positive_selection_results_ready )
    // .set{ all_samples_results }

    PLOTMUTDENSITYQC(all_mutdensities, panel, groups_definition, group_name)
    // mutation density per gene cohort-level
    // mutation density per gene & sample
    //      synonymous
    //      non-protein-affecting
    // pending:
    //      protein-affecting
    //          truncating
    //          missense

    APPLYOMEGAQC(all_omegas, PLOTMUTDENSITYQC.out.compiled_flagged.collect())

    emit:
    mutdensity_plots    = PLOTMUTDENSITYQC.out.plots
    flagged_omegas      = APPLYOMEGAQC.out.all_omegas_annotated

}
