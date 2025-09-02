
include { PLOT_MUTDENSITY_QC                   as PLOTMUTDENSITYQC                    } from '../../../modules/local/plot/qc/mutation_densities/main'



workflow PLOTTING_QC {

    take:
    // positive_selection_results_ready
    all_mutdensities
    // all_samples_depth
    samples
    // all_groups
    panel
    // full_panel_rich
    // seqinfo_df
    // domain_df
    // exons_depths_df


    main:

    // pdb_tool_df   = params.annotations3d
    //                         ? Channel.fromPath( "${params.annotations3d}/pdb_tool_df.tsv", checkIfExists: true).first()
                            // : Channel.empty()


    // plotting only for the entire cohort group
    // Channel.of([ [ id: "all_samples" ] ])
    // .join( positive_selection_results_ready )
    // .set{ all_samples_results }

    PLOTMUTDENSITYQC(all_mutdensities, panel, samples)
    // mutation density per gene cohort-level
    // mutation density per gene & sample
    //      synonymous
    //      non-protein-affecting
    // pending:
    //      protein-affecting
    //          truncating
    //          missense


    emit:
    mutdensity_plots    = PLOTMUTDENSITYQC.out.plots

}
